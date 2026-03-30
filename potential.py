#!/usr/bin/env python3
"""
Calculate electrostatic potential along a box axis from MD simulations.

Reimplements the functionality of GROMACS 'gmx potential' using MDAnalysis.
Computes charge density, electric field, and electrostatic potential by
slicing the system along a chosen axis and solving the Poisson equation.

By default, solves Poisson's equation in Fourier space, which naturally
respects periodic boundary conditions and avoids artifacts from molecule
splitting at the box boundary. Classical real-space double integration
(as in gmx potential) is available via -classical, with an optional
Sachs et al. correction (-sachs).

Reference: Gurtovenko & Vattulainen, J. Chem. Phys. 130, 215107 (2009).
"""

import argparse
import sys
import warnings
import numpy as np
import MDAnalysis as mda

# Physical constants
ELEMENTARY_CHARGE = 1.602176634e-19      # C
EPSILON_0 = 8.8541878128e-12             # C^2 / (N m^2)
NM_TO_M = 1e-9

# Conversion factor: e/nm^3 -> C/m^3
CHARGE_DENSITY_TO_SI = ELEMENTARY_CHARGE / (NM_TO_M ** 3)


def parse_ndx(filename):
    """Parse a GROMACS index (.ndx) file.

    Returns an ordered list of (group_name, atom_indices) tuples.
    Atom indices are 0-based (converted from GROMACS 1-based).
    """
    groups = []
    current_name = None
    current_indices = []

    with open(filename) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith(';'):
                continue
            if line.startswith('['):
                # Save previous group
                if current_name is not None:
                    groups.append((current_name, np.array(current_indices, dtype=int)))
                current_name = line.strip('[] ')
                current_indices = []
            else:
                # Parse atom indices (1-based in .ndx, convert to 0-based)
                current_indices.extend(int(x) - 1 for x in line.split())

    # Save last group
    if current_name is not None:
        groups.append((current_name, np.array(current_indices, dtype=int)))

    return groups


def generate_default_groups(universe):
    """Generate default index groups from the topology, similar to gmx make_ndx.

    Returns an ordered list of (group_name, atom_indices) tuples (0-based).
    """
    groups = []

    # System: all atoms
    groups.append(("System", np.arange(len(universe.atoms))))

    # Per-residue-name groups
    resnames = []
    seen = set()
    for res in universe.atoms.residues:
        if res.resname not in seen:
            resnames.append(res.resname)
            seen.add(res.resname)

    for rn in resnames:
        sel = universe.select_atoms(f"resname {rn}")
        groups.append((rn, sel.ix.copy()))

    return groups


def prompt_group_selection(groups, prompt_text):
    """Display available groups and prompt the user to select one.

    Returns the (group_name, atom_indices) tuple for the chosen group.
    """
    print(f"\n{prompt_text}")
    print("-" * 50)
    for i, (name, indices) in enumerate(groups):
        print(f"  {i:3d}  {name} ({len(indices)} atoms)")
    print("-" * 50)

    while True:
        try:
            choice = input(f"Select group (0-{len(groups)-1}): ").strip()
            idx = int(choice)
            if 0 <= idx < len(groups):
                name, indices = groups[idx]
                print(f"  -> Selected: {name} ({len(indices)} atoms)")
                return groups[idx]
            else:
                print(f"  Invalid choice. Enter a number between 0 and {len(groups)-1}.")
        except ValueError:
            print(f"  Invalid input. Enter a number.")
        except (EOFError, KeyboardInterrupt):
            print("\nAborted.")
            sys.exit(1)


def estimate_n_slices(box_length_nm, n_atoms):
    """Estimate a reasonable number of slices."""
    n_from_resolution = int(box_length_nm / 0.02)
    n_from_atoms = max(10, n_atoms // 5)
    n_slices = min(n_from_resolution, n_from_atoms)
    n_slices = max(50, min(1000, n_slices))
    return n_slices


def write_xvg(filename, x, datasets, title, xlabel, ylabel, legends,
              xvg_format="xmgrace"):
    """Write data in GROMACS .xvg format."""
    with open(filename, "w") as f:
        f.write(f"# This file was created by potential.py\n")
        f.write(f"# Electrostatic potential calculation using MDAnalysis\n")
        if xvg_format != "none":
            f.write(f'@    title "{title}"\n')
            f.write(f'@    xaxis  label "{xlabel}"\n')
            f.write(f'@    yaxis  label "{ylabel}"\n')
            f.write(f"@TYPE xy\n")
            for i, legend in enumerate(legends):
                f.write(f'@ s{i} legend "{legend}"\n')
        for i in range(len(x)):
            line = f"  {x[i]:.6f}"
            for ds in datasets:
                line += f"  {ds[i]:.6f}"
            f.write(line + "\n")


def integrate_poisson_fourier(charge_density_SI, n_slices, slice_width_m):
    """Solve Poisson equation in Fourier space (inherently periodic).

    Solves d^2 psi / dz^2 = -rho / epsilon_0 in reciprocal space:
        psi(k) = rho(k) / (k^2 * epsilon_0)   for k != 0

    The k=0 component (overall potential offset) is set to zero.
    """
    rho_fft = np.fft.rfft(charge_density_SI)
    k = 2 * np.pi * np.fft.rfftfreq(n_slices, d=slice_width_m)

    pot_fft = np.zeros_like(rho_fft)
    pot_fft[1:] = rho_fft[1:] / (k[1:]**2 * EPSILON_0)

    potential = np.fft.irfft(pot_fft, n=n_slices)

    E_fft = np.zeros_like(rho_fft)
    E_fft[1:] = -1j * k[1:] * pot_fft[1:]
    electric_field = np.fft.irfft(E_fft, n=n_slices)

    return potential, electric_field


def integrate_poisson_classical(charge_density_SI, n_slices, slice_width_m,
                                sachs=False, correct=False):
    """Solve Poisson equation by real-space double integration.

    Classical approach as in gmx potential. Uses the gmx potential sign
    convention: psi = +(1/epsilon_0) * double_integral(rho).

    Boundary conditions: psi(0) = 0, E(0) = 0.

    If correct=True, applies the gmx potential -correct procedure:
    subtract the mean charge density before integration, then subtract
    the mean electric field after the first integration. This removes
    both constant and linear drift terms from the potential.

    If sachs=True, applies the Sachs et al. correction instead.
    """
    rho = charge_density_SI.copy()

    if correct:
        # Step 1: subtract mean charge density (removes constant E drift)
        nonzero = np.abs(rho) > np.finfo(float).tiny
        if nonzero.any():
            rho[nonzero] -= rho[nonzero].mean()

    electric_field = np.zeros(n_slices, dtype=np.float64)
    for i in range(1, n_slices):
        electric_field[i] = electric_field[i - 1] + (
            (rho[i - 1] + rho[i]) / 2.0
            * slice_width_m / EPSILON_0
        )

    if correct:
        # Step 2: subtract mean electric field (removes linear psi drift)
        nonzero_E = np.abs(charge_density_SI) > np.finfo(float).tiny
        if nonzero_E.any():
            electric_field[nonzero_E] -= electric_field[nonzero_E].mean()

    potential = np.zeros(n_slices, dtype=np.float64)
    for i in range(1, n_slices):
        potential[i] = potential[i - 1] - (
            (electric_field[i - 1] + electric_field[i]) / 2.0
            * slice_width_m
        )

    if sachs:
        z_frac = (np.arange(n_slices) + 0.5) / n_slices
        potential = potential - z_frac * potential[-1]

    return potential, electric_field


def detect_water_regions(water_density, z_nm, threshold_frac=0.5):
    """Detect bulk water regions from water density profile.

    Returns a boolean mask over slices where water density exceeds
    threshold_frac * max(water_density). Contiguous regions shorter
    than 5 slices are discarded.
    """
    if water_density is None or water_density.max() == 0:
        return None
    threshold = threshold_frac * water_density.max()
    mask = water_density > threshold

    # Remove short contiguous regions (likely not bulk water)
    labeled = np.zeros_like(mask, dtype=int)
    label = 0
    for i in range(len(mask)):
        if mask[i]:
            if i == 0 or not mask[i - 1]:
                label += 1
            labeled[i] = label

    for lbl in range(1, label + 1):
        region = labeled == lbl
        if region.sum() < 5:
            mask[region] = False

    return mask


def find_contiguous_regions(mask, z_nm):
    """Find contiguous True regions in a boolean mask.

    Returns list of (start_idx, end_idx) tuples.
    """
    regions = []
    in_region = False
    for i in range(len(mask)):
        if mask[i] and not in_region:
            start = i
            in_region = True
        elif not mask[i] and in_region:
            regions.append((start, i - 1))
            in_region = False
    if in_region:
        regions.append((start, len(mask) - 1))
    return regions


def measure_water_slope(potential, z_nm, water_mask):
    """Measure the slope of the potential in bulk water regions.

    Fits each contiguous water region separately and returns the
    average slope, since different regions may be at different
    absolute potential levels (the reaction potential is periodic).

    Returns (slope_V_per_nm, slopes_per_region, r_squared_per_region).
    """
    if water_mask is None or water_mask.sum() < 3:
        return None, None, None

    regions = find_contiguous_regions(water_mask, z_nm)
    if not regions:
        return None, None, None

    slopes = []
    r_squareds = []
    for start, end in regions:
        if end - start < 2:
            continue
        z_reg = z_nm[start:end + 1]
        pot_reg = potential[start:end + 1]
        coeffs = np.polyfit(z_reg, pot_reg, 1)
        slopes.append(coeffs[0])

        pred = np.polyval(coeffs, z_reg)
        ss_res = np.sum((pot_reg - pred) ** 2)
        ss_tot = np.sum((pot_reg - pot_reg.mean()) ** 2)
        r_squareds.append(1 - ss_res / ss_tot if ss_tot > 0 else 0.0)

    if not slopes:
        return None, None, None

    avg_slope = np.mean(slopes)
    return avg_slope, slopes, r_squareds


def integrate_poisson_2d_fourier(charge_density_2d_SI, n1, n2, d1_m, d2_m):
    """Solve 2D Poisson equation in Fourier space.

    Solves d^2 psi/dx1^2 + d^2 psi/dx2^2 = -rho / epsilon_0
    in reciprocal space:
        psi(k1,k2) = rho(k1,k2) / ((k1^2 + k2^2) * epsilon_0)
    for (k1,k2) != (0,0).

    Parameters
    ----------
    charge_density_2d_SI : ndarray, shape (n1, n2)
        Charge density in C/m^3, averaged over the third axis.
    n1, n2 : int
        Number of bins along each axis.
    d1_m, d2_m : float
        Bin width along each axis in meters.

    Returns
    -------
    potential : ndarray, shape (n1, n2)
        Electrostatic potential in V.
    """
    rho_fft = np.fft.rfft2(charge_density_2d_SI)

    k1 = 2 * np.pi * np.fft.fftfreq(n1, d=d1_m)
    k2 = 2 * np.pi * np.fft.rfftfreq(n2, d=d2_m)

    K1, K2 = np.meshgrid(k1, k2, indexing='ij')
    k_sq = K1**2 + K2**2

    pot_fft = np.zeros_like(rho_fft)
    nonzero = k_sq > 0
    pot_fft[nonzero] = rho_fft[nonzero] / (k_sq[nonzero] * EPSILON_0)

    potential = np.fft.irfft2(pot_fft, s=(n1, n2))
    return potential


def write_2d_map(filename, x1, x2, data, title, label1, label2, zlabel):
    """Write 2D data in gnuplot-compatible format.

    Format: three columns (axis1, axis2, value) with blank lines
    between rows (gnuplot pm3d / splot format). Also readable by
    numpy.loadtxt after stripping blanks.
    """
    with open(filename, "w") as f:
        f.write(f"# {title}\n")
        f.write(f"# Columns: {label1}  {label2}  {zlabel}\n")
        f.write(f"# Grid: {len(x1)} x {len(x2)}\n")
        for i in range(len(x1)):
            for j in range(len(x2)):
                f.write(f"  {x1[i]:.6f}  {x2[j]:.6f}  {data[i, j]:.6e}\n")
            f.write("\n")


def estimate_n_bins_2d(box_length_nm, n_atoms):
    """Estimate number of bins per axis for 2D maps.

    More conservative than 1D since each bin gets fewer atoms
    (n_atoms / n_bins^2 instead of n_atoms / n_bins).
    """
    n_from_resolution = int(box_length_nm / 0.05)
    n_from_atoms = max(10, int(np.sqrt(n_atoms / 5)))
    n_bins = min(n_from_resolution, n_from_atoms)
    n_bins = max(20, min(500, n_bins))
    return n_bins


def compute_potential_2d(tpr_file, traj_file, plane="XZ",
                         output_potential="potential_2d.dat",
                         output_charge="charge_2d.dat",
                         output_total="potential_total_2d.dat",
                         n_slices=None, begin=None, end=None, dt=None,
                         center=False, group=None, center_group=None,
                         ndx_file=None, efield=None, xvg_format="xmgrace"):
    """
    Compute 2D electrostatic potential map from MD trajectory.

    Bins charges into a 2D grid on the specified plane, averages over
    the third axis, and solves the 2D Poisson equation in Fourier space.

    Parameters
    ----------
    tpr_file : str
        Path to GROMACS .tpr file.
    traj_file : str
        Path to trajectory file.
    plane : str
        Two-letter plane specification: "XZ", "XY", or "YZ".
        First letter = axis 1 (rows), second = axis 2 (columns).
        Charges are averaged over the remaining axis.
    output_potential : str
        Output file for 2D potential map.
    output_charge : str
        Output file for 2D charge density map.
    output_total : str
        Output file for 2D total potential (reaction + applied field ramp).
        Only written when efield is set.
    n_slices : int or None
        Number of bins per axis. If None, auto-estimated.
    begin, end, dt : float or None
        Time filtering options (ps).
    center : bool
        Center coordinates w.r.t. center of mass each frame.
    group : str or None
        MDAnalysis atom selection.
    center_group : str or None
        MDAnalysis selection for centering.
    ndx_file : str or None
        GROMACS index file for interactive group selection.
    efield : list of 3 floats, or None
        Applied electric field [Ex, Ey, Ez] in V/nm. When set, computes
        the total potential by adding the linear ramp from the in-plane
        field components to the reaction potential.
    xvg_format : str
        Unused for 2D output, kept for interface consistency.
    """
    axis_map = {"X": 0, "Y": 1, "Z": 2}
    plane = plane.upper()
    if len(plane) != 2 or plane[0] not in axis_map or plane[1] not in axis_map:
        raise ValueError(f"Invalid plane '{plane}'. Must be two of X, Y, Z "
                         f"(e.g., XZ, XY, YZ).")
    if plane[0] == plane[1]:
        raise ValueError(f"Plane axes must be different, got '{plane}'.")

    ax1_idx = axis_map[plane[0]]
    ax2_idx = axis_map[plane[1]]
    avg_axes = [i for i in range(3) if i not in (ax1_idx, ax2_idx)]
    avg_ax_idx = avg_axes[0]
    avg_ax_name = {0: "X", 1: "Y", 2: "Z"}[avg_ax_idx]

    print(f"2D potential map on {plane[0]}{plane[1]} plane "
          f"(averaging over {avg_ax_name})")

    # Load universe
    print(f"Loading topology from {tpr_file}...")
    u = mda.Universe(tpr_file, traj_file)

    if not hasattr(u.atoms, "charges"):
        raise RuntimeError(
            "No charges found in the topology file. "
            "Ensure you are using a .tpr file that contains charge information."
        )

    # --- Group selection ---
    if ndx_file is not None:
        groups = parse_ndx(ndx_file)
        if not groups:
            raise RuntimeError(f"No groups found in {ndx_file}")
        calc_name, calc_indices = prompt_group_selection(
            groups, "Select group for potential calculation:")
        atoms = u.atoms[calc_indices]
        if center:
            center_name, center_indices = prompt_group_selection(
                groups, "Select group for centering:")
            center_atoms = u.atoms[center_indices]
        else:
            center_atoms = atoms
    else:
        if group is None:
            group = "all"
        atoms = u.select_atoms(group)
        if center_group is None:
            center_atoms = atoms
        else:
            center_atoms = u.select_atoms(center_group)

    print(f"Selected {len(atoms)} atoms for analysis.")
    if center:
        print(f"Centering group: {len(center_atoms)} atoms.")

    charges = atoms.charges

    # Estimate grid size
    u.trajectory[0]
    box = u.dimensions[:3]
    box1_nm = box[ax1_idx] / 10.0
    box2_nm = box[ax2_idx] / 10.0

    if n_slices is None:
        n1 = estimate_n_bins_2d(box1_nm, len(atoms))
        n2 = estimate_n_bins_2d(box2_nm, len(atoms))
        print(f"Auto-selected {n1} x {n2} bins "
              f"(box: {box1_nm:.2f} x {box2_nm:.2f} nm).")
    else:
        n1 = n_slices
        n2 = n_slices
        print(f"Using {n1} x {n2} bins.")

    print(f"Integration method: 2D Fourier")

    if efield is not None:
        e1 = efield[ax1_idx]
        e2 = efield[ax2_idx]
        e_avg = efield[avg_ax_idx]
        print(f"Applied electric field: E = [{efield[0]}, {efield[1]}, {efield[2]}] V/nm")
        print(f"  In-plane components: E_{plane[0]} = {e1} V/nm, "
              f"E_{plane[1]} = {e2} V/nm")
        if e_avg != 0:
            print(f"  Note: E_{avg_ax_name} = {e_avg} V/nm component is along "
                  f"the averaging axis and does not appear in the 2D map.")

    # Accumulate charge density over frames
    charge_density_sum = np.zeros((n1, n2), dtype=np.float64)
    n_frames = 0
    box1_sum = 0.0
    box2_sum = 0.0
    box_avg_sum = 0.0

    print("Processing trajectory...")
    for ts in u.trajectory:
        t = ts.time
        if begin is not None and t < begin:
            continue
        if end is not None and t > end:
            break
        if dt is not None and n_frames > 0:
            if abs(t % dt) > 0.01 and abs(t % dt - dt) > 0.01:
                continue

        box = ts.dimensions[:3]
        box1 = box[ax1_idx]
        box2 = box[ax2_idx]
        box_avg = box[avg_ax_idx]

        pos1 = atoms.positions[:, ax1_idx]
        pos2 = atoms.positions[:, ax2_idx]

        if center:
            com = center_atoms.center_of_mass()
            pos1 = pos1 - com[ax1_idx] + box1 / 2.0
            pos2 = pos2 - com[ax2_idx] + box2 / 2.0

        pos1 = pos1 % box1
        pos2 = pos2 % box2

        d1 = box1 / n1
        d2 = box2 / n2
        idx1 = np.floor(pos1 / d1).astype(int)
        idx2 = np.floor(pos2 / d2).astype(int)
        idx1 = np.clip(idx1, 0, n1 - 1)
        idx2 = np.clip(idx2, 0, n2 - 1)

        # Bin volume: d1 * d2 * box_avg (full extent of averaging axis)
        bin_volume = d1 * d2 * box_avg

        frame_charge = np.zeros((n1, n2), dtype=np.float64)
        linear_idx = idx1 * n2 + idx2
        np.add.at(frame_charge.ravel(), linear_idx, charges)

        charge_density_sum += frame_charge / bin_volume

        box1_sum += box1
        box2_sum += box2
        box_avg_sum += box_avg
        n_frames += 1

    if n_frames == 0:
        raise RuntimeError("No frames were processed. Check -b/-e/-dt options.")

    print(f"Processed {n_frames} frames.")

    charge_density = charge_density_sum / n_frames

    avg_box1 = box1_sum / n_frames
    avg_box2 = box2_sum / n_frames
    avg_box1_nm = avg_box1 / 10.0
    avg_box2_nm = avg_box2 / 10.0

    # Convert: charge_density is in e/A^3, convert to e/nm^3 then to C/m^3
    charge_density_e_nm3 = charge_density * 1e3  # e/A^3 -> e/nm^3
    charge_density_SI = charge_density_e_nm3 * CHARGE_DENSITY_TO_SI

    d1_nm = avg_box1_nm / n1
    d2_nm = avg_box2_nm / n2
    d1_m = d1_nm * NM_TO_M
    d2_m = d2_nm * NM_TO_M

    # Coordinates at bin centers
    x1_nm = (np.arange(n1) + 0.5) * d1_nm
    x2_nm = (np.arange(n2) + 0.5) * d2_nm

    # Solve 2D Poisson
    potential = integrate_poisson_2d_fourier(
        charge_density_SI, n1, n2, d1_m, d2_m)

    # Write outputs
    write_2d_map(output_charge, x1_nm, x2_nm, charge_density_e_nm3,
                 f"2D charge density on {plane[0]}{plane[1]} plane",
                 f"{plane[0]} (nm)", f"{plane[1]} (nm)",
                 "Charge density (e/nm^3)")
    print(f"Wrote 2D charge density to {output_charge}")

    write_2d_map(output_potential, x1_nm, x2_nm, potential,
                 f"2D electrostatic potential on {plane[0]}{plane[1]} plane",
                 f"{plane[0]} (nm)", f"{plane[1]} (nm)",
                 "Potential (V)")
    print(f"Wrote 2D potential to {output_potential}")

    # --- Applied electric field: total potential ---
    if efield is not None:
        e1 = efield[ax1_idx]
        e2 = efield[ax2_idx]

        # Total potential = reaction potential - E1*x1 - E2*x2
        X1, X2 = np.meshgrid(x1_nm, x2_nm, indexing='ij')
        total_potential = potential - e1 * X1 - e2 * X2

        write_2d_map(output_total, x1_nm, x2_nm, total_potential,
                     f"2D total potential (reaction + applied field) "
                     f"on {plane[0]}{plane[1]} plane",
                     f"{plane[0]} (nm)", f"{plane[1]} (nm)",
                     "Potential (V)")
        print(f"Wrote 2D total potential to {output_total}")

        # Report voltages from in-plane field components
        print(f"\n--- Applied electric field analysis ---")
        print(f"Applied field: E = [{efield[0]}, {efield[1]}, {efield[2]}] V/nm")
        if e1 != 0:
            v1 = e1 * avg_box1_nm
            print(f"Voltage along {plane[0]}: E_{plane[0]} * L_{plane[0]} = "
                  f"{v1*1000:.1f} mV")
        if e2 != 0:
            v2 = e2 * avg_box2_nm
            print(f"Voltage along {plane[1]}: E_{plane[1]} * L_{plane[1]} = "
                  f"{v2*1000:.1f} mV")
        print(f"Total potential range: [{total_potential.min():.4f}, "
              f"{total_potential.max():.4f}] V")

    # Summary
    print(f"\n--- 2D Map Summary ---")
    print(f"Plane: {plane[0]}{plane[1]} (averaged over {avg_ax_name})")
    print(f"Box dimensions: {avg_box1_nm:.3f} x {avg_box2_nm:.3f} nm")
    print(f"Grid: {n1} x {n2} bins")
    print(f"Bin size: {d1_nm:.4f} x {d2_nm:.4f} nm")
    print(f"Frames analyzed: {n_frames}")
    print(f"Potential range: [{potential.min():.4f}, {potential.max():.4f}] V")
    print(f"Charge density range: [{charge_density_e_nm3.min():.4f}, "
          f"{charge_density_e_nm3.max():.4f}] e/nm^3")


def compute_potential(tpr_file, traj_file, output_potential="potential.xvg",
                      output_charge="charge.xvg", output_field="field.xvg",
                      output_total="potential_total.xvg",
                      axis="Z", n_slices=None, begin=None, end=None, dt=None,
                      center=False, symmetrize=False, group=None,
                      center_group=None, ndx_file=None, classical=False,
                      sachs=False, correct=False, efield=None,
                      xvg_format="xmgrace"):
    """
    Compute electrostatic potential from MD trajectory.

    Parameters
    ----------
    tpr_file : str
        Path to GROMACS .tpr file (provides charges and topology).
    traj_file : str
        Path to trajectory file (.xtc, .trr, .gro, etc.).
    output_potential : str
        Output file for electrostatic potential.
    output_charge : str
        Output file for charge density.
    output_field : str
        Output file for electric field.
    output_total : str
        Output file for total potential (reaction + applied field ramp).
        Only written when efield is set.
    axis : str
        Axis along which to compute the potential ('X', 'Y', or 'Z').
    n_slices : int or None
        Number of slices. If None, an automatic estimate is used.
    begin : float or None
        First time to read (ps).
    end : float or None
        Last time to read (ps).
    dt : float or None
        Only use frames at multiples of dt (ps).
    center : bool
        Center coordinates w.r.t. center of mass of center_group each frame.
    symmetrize : bool
        Symmetrize the profiles around the box center (implies center=True).
    group : str or None
        MDAnalysis atom selection for the atoms to include, or None to use
        index group selection.
    center_group : str or None
        MDAnalysis atom selection for centering, or None.
    ndx_file : str or None
        Path to GROMACS index file (.ndx). If provided, groups are selected
        interactively from this file. Overrides -group and -center-group.
    classical : bool
        Use classical real-space double integration instead of Fourier.
    sachs : bool
        With -classical: apply the Sachs et al. correction.
    correct : bool
        With -classical: apply the gmx potential -correct procedure
        (subtract mean charge density, then mean electric field).
    efield : float or None
        Applied external electric field in V/nm (from GROMACS mdp).
        When set, computes total potential (reaction + linear ramp),
        detects water regions, and reports the applied voltage.
    xvg_format : str
        XVG format: "xmgrace", "xmgr", or "none".
    """
    axis_map = {"X": 0, "Y": 1, "Z": 2}
    if axis.upper() not in axis_map:
        raise ValueError(f"Invalid axis '{axis}'. Must be X, Y, or Z.")
    ax_idx = axis_map[axis.upper()]

    if symmetrize:
        center = True

    # Load universe
    print(f"Loading topology from {tpr_file}...")
    u = mda.Universe(tpr_file, traj_file)

    # Check that charges are available
    if not hasattr(u.atoms, "charges"):
        raise RuntimeError(
            "No charges found in the topology file. "
            "Ensure you are using a .tpr file that contains charge information."
        )

    # --- Group selection ---
    if ndx_file is not None:
        # Load groups from index file
        groups = parse_ndx(ndx_file)
        if not groups:
            raise RuntimeError(f"No groups found in {ndx_file}")

        # Select group for potential calculation
        calc_name, calc_indices = prompt_group_selection(
            groups, "Select group for potential calculation:")
        atoms = u.atoms[calc_indices]

        # Select group for centering
        if center or symmetrize:
            center_name, center_indices = prompt_group_selection(
                groups, "Select group for centering:")
            center_atoms = u.atoms[center_indices]
        else:
            center_atoms = atoms
    else:
        # Use MDAnalysis selection strings
        if group is None:
            group = "all"
        atoms = u.select_atoms(group)

        if center_group is None:
            center_atoms = atoms
        else:
            center_atoms = u.select_atoms(center_group)

    print(f"Selected {len(atoms)} atoms for analysis.")
    if center or symmetrize:
        print(f"Centering group: {len(center_atoms)} atoms.")

    charges = atoms.charges  # in elementary charge units

    # Estimate number of slices if not provided
    if n_slices is None:
        u.trajectory[0]
        box_length = u.dimensions[ax_idx]
        box_length_nm = box_length / 10.0
        n_slices = estimate_n_slices(box_length_nm, len(atoms))
        print(f"Auto-selected {n_slices} slices "
              f"(box length along {axis.upper()}: {box_length_nm:.2f} nm).")
    else:
        print(f"Using {n_slices} slices.")

    method = "classical" if classical else "Fourier"
    if classical and sachs:
        method += " + Sachs correction"
    if classical and correct:
        method += " + correct"
    print(f"Integration method: {method}")

    if efield is not None:
        print(f"Applied electric field: {efield} V/nm")

    # Select water atoms for bulk water detection (used with -efield)
    water_atoms = None
    if efield is not None:
        for water_sel in ["resname SOL", "resname TIP3", "resname HOH",
                          "resname WAT", "resname SPC"]:
            try:
                sel = u.select_atoms(water_sel)
                if len(sel) > 0:
                    water_atoms = sel
                    print(f"Water detection: {len(water_atoms)} atoms "
                          f"(resname {water_sel.split()[-1]})")
                    break
            except Exception:
                continue
        if water_atoms is None:
            print("Warning: could not detect water atoms for bulk region "
                  "identification. Voltage estimation from slope will be skipped.")

    # Accumulate charge density histogram over frames
    charge_density_sum = np.zeros(n_slices, dtype=np.float64)
    water_density_sum = np.zeros(n_slices, dtype=np.float64) if water_atoms is not None else None
    n_frames = 0
    box_length_sum = 0.0
    box_area_sum = 0.0

    cross_axes = [i for i in range(3) if i != ax_idx]

    print("Processing trajectory...")
    for ts in u.trajectory:
        t = ts.time
        if begin is not None and t < begin:
            continue
        if end is not None and t > end:
            break
        if dt is not None and n_frames > 0:
            if abs(t % dt) > 0.01 and abs(t % dt - dt) > 0.01:
                continue

        box = ts.dimensions[:3]
        box_length = box[ax_idx]
        box_area = box[cross_axes[0]] * box[cross_axes[1]]

        positions = atoms.positions[:, ax_idx]

        if center:
            com = center_atoms.center_of_mass()[ax_idx]
            positions = positions - com + box_length / 2.0

        positions = positions % box_length

        slice_width = box_length / n_slices
        slice_indices = np.floor(positions / slice_width).astype(int)
        slice_indices = np.clip(slice_indices, 0, n_slices - 1)

        frame_charge = np.zeros(n_slices, dtype=np.float64)
        np.add.at(frame_charge, slice_indices, charges)

        slice_volume = slice_width * box_area
        frame_charge_density = frame_charge / slice_volume

        charge_density_sum += frame_charge_density

        # Track water density for bulk region detection
        if water_atoms is not None:
            water_pos = water_atoms.positions[:, ax_idx]
            if center:
                water_pos = water_pos - com + box_length / 2.0
            water_pos = water_pos % box_length
            water_slice_idx = np.floor(water_pos / slice_width).astype(int)
            water_slice_idx = np.clip(water_slice_idx, 0, n_slices - 1)
            water_count = np.zeros(n_slices, dtype=np.float64)
            np.add.at(water_count, water_slice_idx, 1.0)
            water_density_sum += water_count / slice_volume

        box_length_sum += box_length
        box_area_sum += box_area
        n_frames += 1

    if n_frames == 0:
        raise RuntimeError("No frames were processed. Check -b/-e/-dt options.")

    print(f"Processed {n_frames} frames.")

    charge_density = charge_density_sum / n_frames
    avg_box_length = box_length_sum / n_frames
    avg_box_length_nm = avg_box_length / 10.0

    charge_density_e_nm3 = charge_density * 1e3
    charge_density_SI = charge_density_e_nm3 * CHARGE_DENSITY_TO_SI

    slice_width_nm = avg_box_length_nm / n_slices
    slice_width_m = slice_width_nm * NM_TO_M

    if symmetrize:
        charge_density_e_nm3 = (charge_density_e_nm3
                                + charge_density_e_nm3[::-1]) / 2.0
        charge_density_SI = (charge_density_SI
                             + charge_density_SI[::-1]) / 2.0

    z_nm = (np.arange(n_slices) + 0.5) * slice_width_nm

    if classical:
        potential, electric_field = integrate_poisson_classical(
            charge_density_SI, n_slices, slice_width_m, sachs=sachs,
            correct=correct)
    else:
        potential, electric_field = integrate_poisson_fourier(
            charge_density_SI, n_slices, slice_width_m)

    electric_field_V_nm = electric_field * NM_TO_M

    if symmetrize:
        electric_field_V_nm = (electric_field_V_nm
                               - electric_field_V_nm[::-1]) / 2.0
        potential = (potential + potential[::-1]) / 2.0

    legends = ["System"]

    write_xvg(output_charge, z_nm, [charge_density_e_nm3],
              "Charge density", f"{axis.upper()} (nm)", "Charge density (e/nm^3)",
              legends, xvg_format)
    print(f"Wrote charge density to {output_charge}")

    write_xvg(output_field, z_nm, [electric_field_V_nm],
              "Electric field", f"{axis.upper()} (nm)", "E (V/nm)",
              legends, xvg_format)
    print(f"Wrote electric field to {output_field}")

    write_xvg(output_potential, z_nm, [potential],
              "Electrostatic potential", f"{axis.upper()} (nm)", "Potential (V)",
              legends, xvg_format)
    print(f"Wrote potential to {output_potential}")

    pot_drop = potential[-1] - potential[0]
    pot_min = np.min(potential)
    pot_max = np.max(potential)
    print(f"\n--- Summary ---")
    print(f"Box length along {axis.upper()}: {avg_box_length_nm:.3f} nm")
    print(f"Number of slices: {n_slices}")
    print(f"Slice width: {slice_width_nm:.4f} nm")
    print(f"Frames analyzed: {n_frames}")
    print(f"Potential range: [{pot_min:.4f}, {pot_max:.4f}] V")
    print(f"Potential drop (last - first): {pot_drop:.4f} V")

    rho_integral = np.sum(charge_density_e_nm3) * slice_width_nm
    if abs(rho_integral) > 0.01:
        warnings.warn(
            f"Integrated charge density along the axis is {rho_integral:.4f} "
            f"e/nm^2, which is not close to zero. This may indicate a net "
            f"charge in the system or insufficient sampling. The potential "
            f"profile may show a linear drift."
        )

    # --- Applied electric field analysis ---
    if efield is not None:
        voltage_exact = efield * avg_box_length_nm
        print(f"\n--- Applied electric field analysis ---")
        print(f"Applied field E_z: {efield} V/nm")
        print(f"Average box length: {avg_box_length_nm:.4f} nm")
        print(f"Applied voltage V = E_z * L_z: {voltage_exact*1000:.1f} mV")

        # Total potential = reaction potential - E * z
        total_potential = potential - efield * z_nm

        # Total electric field = reaction field + applied field
        total_field_V_nm = electric_field_V_nm + efield

        write_xvg(output_total, z_nm, [total_potential],
                  "Total electrostatic potential (reaction + applied field)",
                  f"{axis.upper()} (nm)", "Potential (V)",
                  ["Total"], xvg_format)
        print(f"Wrote total potential to {output_total}")

        # Detect water regions and measure slope
        water_density = water_density_sum / n_frames if water_density_sum is not None else None
        water_mask = detect_water_regions(water_density, z_nm)

        if water_mask is not None and water_mask.sum() >= 3:
            n_water_slices = water_mask.sum()
            z_water_ranges = []
            in_region = False
            for i in range(len(water_mask)):
                if water_mask[i] and not in_region:
                    start = z_nm[i]
                    in_region = True
                elif not water_mask[i] and in_region:
                    z_water_ranges.append((start, z_nm[i - 1]))
                    in_region = False
            if in_region:
                z_water_ranges.append((start, z_nm[-1]))

            print(f"Detected {len(z_water_ranges)} bulk water region(s) "
                  f"({n_water_slices} slices):")
            for r_start, r_end in z_water_ranges:
                print(f"  z = {r_start:.2f} - {r_end:.2f} nm")

            # Measure slope of reaction potential in water
            slope, slopes, r_squareds = measure_water_slope(
                potential, z_nm, water_mask)

            if slope is not None:
                # Report per-region slopes
                print(f"\nReaction potential slope in water (per region):")
                for i, (sl, rsq) in enumerate(zip(slopes, r_squareds)):
                    print(f"  Region {i+1}: slope = {sl:.4f} V/nm "
                          f"(R^2 = {rsq:.4f})")

                voltage_from_slope = slope * avg_box_length_nm
                print(f"Average slope: {slope:.4f} V/nm")
                print(f"Voltage from slope (E_slope * L_z): "
                      f"{voltage_from_slope*1000:.1f} mV")
                print(f"Voltage from E_applied * L_z:       "
                      f"{voltage_exact*1000:.1f} mV")
                recovery = slope / efield * 100 if efield != 0 else 0
                print(f"Recovery: {recovery:.1f}%")

                # Also measure slope of reaction E field in water
                E_water = electric_field_V_nm[water_mask]
                avg_E_water = E_water.mean()
                print(f"\nAverage reaction E field in water: "
                      f"{avg_E_water:.4f} V/nm")
                print(f"Expected (= -E_applied):           "
                      f"{-efield:.4f} V/nm")

                # Total E field in water (should be ~0 for perfect conductor)
                avg_total_E_water = total_field_V_nm[water_mask].mean()
                print(f"Average total E field in water:     "
                      f"{avg_total_E_water:.4f} V/nm (ideal: 0)")
        else:
            print("Could not detect bulk water regions for slope analysis.")


def main():
    parser = argparse.ArgumentParser(
        description="Calculate electrostatic potential along a box axis "
                    "from MD simulations (analogous to gmx potential).",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""\
Examples:
  %(prog)s -s topol.tpr -f traj.xtc -center
  %(prog)s -s topol.tpr -f traj.xtc -n index.ndx -center
  %(prog)s -s topol.tpr -f traj.xtc -d Z -sl 500 -center -group "not resname SOL"
  %(prog)s -s topol.tpr -f traj.xtc -symm -center-group "resname POPC"
  %(prog)s -s topol.tpr -f traj.xtc -classical -sachs -sl 300
"""
    )

    # Input files
    parser.add_argument("-s", dest="tpr", required=True,
                        help="Input topology file (.tpr)")
    parser.add_argument("-f", dest="traj", required=True,
                        help="Input trajectory file (.xtc, .trr, .gro, ...)")
    parser.add_argument("-n", dest="ndx", default=None,
                        help="GROMACS index file (.ndx) for group selection. "
                             "When provided, groups are chosen interactively "
                             "(overrides -group and -center-group).")

    # Output files
    parser.add_argument("-o", dest="output_potential", default="potential.xvg",
                        help="Output: electrostatic potential (default: potential.xvg)")
    parser.add_argument("-oc", dest="output_charge", default="charge.xvg",
                        help="Output: charge density (default: charge.xvg)")
    parser.add_argument("-of", dest="output_field", default="field.xvg",
                        help="Output: electric field (default: field.xvg)")
    parser.add_argument("-ot", dest="output_total", default="potential_total.xvg",
                        help="Output: total potential with applied field "
                             "(default: potential_total.xvg, requires -efield)")

    # Analysis options
    parser.add_argument("-d", dest="axis", default="Z", choices=["X", "Y", "Z"],
                        help="Axis along which to compute the potential (default: Z)")
    parser.add_argument("-sl", dest="n_slices", type=int, default=None,
                        help="Number of slices (default: auto-estimated)")
    parser.add_argument("-b", dest="begin", type=float, default=None,
                        help="First frame time to read (ps)")
    parser.add_argument("-e", dest="end", type=float, default=None,
                        help="Last frame time to read (ps)")
    parser.add_argument("-dt", dest="dt", type=float, default=None,
                        help="Only use frames at multiples of dt (ps)")

    # Centering and symmetrization
    parser.add_argument("-center", action="store_true",
                        help="Center coordinates w.r.t. center of mass of "
                             "center-group each frame")
    parser.add_argument("-symm", action="store_true",
                        help="Symmetrize profiles around center (implies -center)")
    parser.add_argument("-group", dest="group", default=None,
                        help="MDAnalysis selection for atoms to include "
                             '(default: "all"). Ignored if -n is provided.')
    parser.add_argument("-center-group", dest="center_group", default=None,
                        help="MDAnalysis selection for centering group "
                             "(default: same as -group). Ignored if -n is provided.")

    # Integration method
    parser.add_argument("-classical", action="store_true",
                        help="Use classical real-space double integration "
                             "instead of Fourier (reproduces gmx potential behavior)")
    parser.add_argument("-sachs", action="store_true",
                        help="With -classical: apply Sachs et al. correction "
                             "(enforce equal potential on both box sides)")
    parser.add_argument("-correct", action="store_true",
                        help="With -classical: apply the gmx potential -correct "
                             "procedure (subtract mean charge density and mean "
                             "electric field to remove linear drift)")

    # 2D potential maps
    parser.add_argument("-2Dmap", dest="map2d", default=None,
                        metavar="PLANE",
                        help="Compute 2D potential map on the given plane "
                             "(e.g., XZ, XY, YZ). Charges are averaged over "
                             "the third axis. Uses 2D Fourier solver. "
                             "Output: potential_2d.dat, charge_2d.dat.")
    parser.add_argument("-2Defield", dest="efield_2d", nargs=3, type=float,
                        default=None, metavar=("EX", "EY", "EZ"),
                        help="Applied electric field for 2D maps: three "
                             "components Ex Ey Ez in V/nm (from GROMACS mdp "
                             "electric-field-x/y/z). The in-plane components "
                             "are added as a linear ramp to produce the total "
                             "potential (potential_total_2d.dat). Requires "
                             "-2Dmap.")

    # Applied electric field
    parser.add_argument("-efield", dest="efield", type=float, default=None,
                        help="Applied external electric field in V/nm (from "
                             "GROMACS mdp electric-field-z). Computes total "
                             "potential (reaction + applied ramp), detects "
                             "bulk water regions, and reports the applied "
                             "voltage V = E * L_z.")

    # Output format
    parser.add_argument("-xvg", dest="xvg_format", default="xmgrace",
                        choices=["xmgrace", "xmgr", "none"],
                        help="XVG output format (default: xmgrace)")

    args = parser.parse_args()

    if args.sachs and not args.classical:
        print("Note: -sachs has no effect without -classical (Fourier method "
              "is inherently periodic).")

    if args.correct and not args.classical:
        print("Note: -correct has no effect without -classical (Fourier method "
              "does not need drift correction).")

    if args.correct and args.sachs:
        print("Note: -correct and -sachs both remove linear drift. "
              "Using both simultaneously is redundant.")

    if args.ndx and (args.group or args.center_group):
        print("Note: -n overrides -group and -center-group. "
              "Groups will be selected interactively from the index file.")

    if args.efield_2d is not None and args.map2d is None:
        print("Note: -2Defield requires -2Dmap. Ignoring -2Defield.")

    if args.map2d is not None:
        compute_potential_2d(
            tpr_file=args.tpr,
            traj_file=args.traj,
            plane=args.map2d,
            output_potential="potential_2d.dat",
            output_charge="charge_2d.dat",
            output_total="potential_total_2d.dat",
            n_slices=args.n_slices,
            begin=args.begin,
            end=args.end,
            dt=args.dt,
            center=args.center,
            group=args.group,
            center_group=args.center_group,
            ndx_file=args.ndx,
            efield=args.efield_2d,
            xvg_format=args.xvg_format,
        )
        return

    compute_potential(
        tpr_file=args.tpr,
        traj_file=args.traj,
        output_potential=args.output_potential,
        output_charge=args.output_charge,
        output_field=args.output_field,
        output_total=args.output_total,
        axis=args.axis,
        n_slices=args.n_slices,
        begin=args.begin,
        end=args.end,
        dt=args.dt,
        center=args.center,
        symmetrize=args.symm,
        group=args.group,
        center_group=args.center_group,
        ndx_file=args.ndx,
        classical=args.classical,
        sachs=args.sachs,
        correct=args.correct,
        efield=args.efield,
        xvg_format=args.xvg_format,
    )


if __name__ == "__main__":
    main()
