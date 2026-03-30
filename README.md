# potential.py

A Python reimplementation of GROMACS `gmx potential` using MDAnalysis. Calculates the electrostatic potential along a simulation box axis by dividing the system into slices, computing the charge density, and solving the Poisson equation.

By default, uses **Fourier-space integration**, which naturally respects periodic boundary conditions and eliminates the slice-count-dependent artifacts that plague the classical real-space double integration used by `gmx potential`.

Based on the methodology described in:
> Gurtovenko & Vattulainen, *J. Chem. Phys.* **130**, 215107 (2009) — [doi:10.1063/1.3148885](https://doi.org/10.1063/1.3148885)

## Requirements

- Python 3.8+
- NumPy
- MDAnalysis (`pip install MDAnalysis`)

Install dependencies:

```bash
pip install numpy MDAnalysis
```

## Quick start

```bash
python potential.py -s topol.tpr -f traj.xtc -center
```

This reads charges from `topol.tpr`, processes every frame of `traj.xtc`, slices the box along the Z axis, and writes three output files:

| File | Contents | Units |
|------|----------|-------|
| `potential.xvg` | Electrostatic potential | V |
| `charge.xvg` | Charge density | e/nm^3 |
| `field.xvg` | Electric field | V/nm |

## Command-line options

### Input / output

| Flag | Default | Description |
|------|---------|-------------|
| `-s` | *(required)* | GROMACS `.tpr` file (provides atomic charges) |
| `-f` | *(required)* | Trajectory file (`.xtc`, `.trr`, `.gro`, `.pdb`, ...) |
| `-o` | `potential.xvg` | Output electrostatic potential (reaction potential from charges) |
| `-oc` | `charge.xvg` | Output charge density |
| `-of` | `field.xvg` | Output electric field |
| `-ot` | `potential_total.xvg` | Output total potential including applied field ramp (requires `-efield`) |

### Analysis

| Flag | Default | Description |
|------|---------|-------------|
| `-d` | `Z` | Axis to slice along (`X`, `Y`, or `Z`) |
| `-sl` | auto | Number of slices (see [Choosing the number of slices](#choosing-the-number-of-slices)) |
| `-b` | -- | First frame time to read (ps) |
| `-e` | -- | Last frame time to read (ps) |
| `-dt` | -- | Only use frames at multiples of this interval (ps) |

### Centering and symmetrization

| Flag | Description |
|------|-------------|
| `-center` | Center coordinates each frame w.r.t. center of mass of the centering group. Eliminates artifacts from bilayer drift along the box axis. |
| `-symm` | Symmetrize profiles around the box center. Implies `-center`. Useful for symmetric bilayers to improve statistics. |
| `-group` | MDAnalysis atom selection string for atoms to include in the calculation (default: `"all"`). |
| `-center-group` | MDAnalysis atom selection string for the centering group (default: same as `-group`). |

### Integration method

| Flag | Description |
|------|-------------|
| `-classical` | Use classical real-space double integration instead of Fourier. Reproduces `gmx potential` behavior, including its artifacts. |
| `-sachs` | With `-classical`: apply the Sachs et al. correction, enforcing equal potential on both sides of the box. Has no effect without `-classical` (Fourier is inherently periodic). |
| `-correct` | With `-classical`: apply the `gmx potential` `-correct` procedure (subtract mean charge density and mean electric field). Removes linear drift. Has no effect without `-classical`. |

### 2D potential maps

| Flag | Default | Description |
|------|---------|-------------|
| `-2Dmap` | -- | Compute a 2D potential map on the given plane. Takes a two-letter plane specification: `XZ`, `XY`, or `YZ`. Charges are averaged over the third axis and the 2D Poisson equation is solved in Fourier space. Output: `potential_2d.dat` and `charge_2d.dat`. |
| `-2Defield` | -- | Applied electric field for 2D maps: three components `Ex Ey Ez` in V/nm (from GROMACS `electric-field-x/y/z` mdp parameters). The in-plane field components are added as a linear ramp to produce the total potential (`potential_total_2d.dat`). Requires `-2Dmap`. |

When `-2Dmap` is used, the tool switches to 2D mode. The `-sl` flag sets the number of bins per axis (same for both). If not specified, the grid size is auto-estimated (more conservatively than 1D, typically 20-500 bins per axis, since each bin receives fewer atoms). The `-center` flag centers along both in-plane axes. The output uses a gnuplot-compatible 3-column format (axis1, axis2, value) with blank lines between rows.

When `-2Defield Ex Ey Ez` is provided, the tool computes the total potential by adding the linear ramp from the in-plane field components: Phi_total(x1, x2) = Phi_reaction(x1, x2) - E1 * x1 - E2 * x2, where E1 and E2 are the field components along the two in-plane axes. The field component along the averaging axis does not appear in the 2D map (the tool prints a note if it is nonzero). The applied voltages along each in-plane axis are reported as V = E * L.

The 1D-specific flags (`-classical`, `-sachs`, `-correct`, `-symm`, `-efield`, `-o`, `-oc`, `-of`, `-ot`) have no effect in 2D mode — the 2D solver always uses Fourier integration.

### Applied electric field

| Flag | Default | Description |
|------|---------|-------------|
| `-efield` | -- | Applied external electric field in V/nm (from GROMACS `electric-field-z` mdp parameter). Computes total potential (reaction + linear ramp), automatically detects bulk water regions, and reports the applied voltage V = E * L_z along with a slope-based voltage estimate. |

### Output format

| Flag | Default | Description |
|------|---------|-------------|
| `-xvg` | `xmgrace` | XVG format: `xmgrace`, `xmgr`, or `none` (plain columns) |

## Usage examples

### Basic usage (Fourier, recommended)

```bash
python potential.py -s topol.tpr -f traj.xtc -center
```

### Classical integration (gmx potential behavior)

```bash
python potential.py -s topol.tpr -f traj.xtc -center -classical
```

### Classical with -correct (drift removal)

```bash
python potential.py -s topol.tpr -f traj.xtc -center -classical -correct
```

### Classical with Sachs correction (symmetric systems)

```bash
python potential.py -s topol.tpr -f traj.xtc -center -classical -sachs
```

### Lipid bilayer with centering and symmetrization

```bash
python potential.py -s topol.tpr -f traj.xtc -center -symm \
    -center-group "resname POPC POPE"
```

### Compute potential for lipids only (excluding solvent)

```bash
python potential.py -s topol.tpr -f traj.xtc -group "not resname SOL NA CL" \
    -center -center-group "resname POPC"
```

### 2D potential map on the XZ plane

```bash
python potential.py -s topol.tpr -f traj.xtc -center -2Dmap XZ
```

This bins charges into a 2D grid on the XZ plane (averaging over Y), solves the 2D Poisson equation in Fourier space, and writes `potential_2d.dat` and `charge_2d.dat`. The output is in gnuplot `splot` format (3 columns: X, Z, value). To visualize with gnuplot:

```gnuplot
set pm3d map
splot 'potential_2d.dat' using 1:2:3 with pm3d title "Potential (V)"
```

### 2D map with custom resolution

```bash
python potential.py -s topol.tpr -f traj.xtc -center -2Dmap XZ -sl 100
```

Uses 100 x 100 bins. Coarser grids converge faster (fewer frames needed); finer grids require longer trajectories.

### 2D map with applied electric field

```bash
python potential.py -s topol.tpr -f traj.xtc -center -2Dmap XZ \
    -2Defield 0 0 E_z
```

Computes the 2D reaction potential on the XZ plane and adds the linear ramp from the applied field (E_z) to produce the total potential in `potential_total_2d.dat`. The total potential shows the non-periodic voltage drop across the membrane, resolved laterally — useful for visualizing how the potential landscape differs inside a channel pore vs the surrounding lipid.

### Simulation with applied electric field

```bash
python potential.py -s topol.tpr -f traj.xtc -center -efield E_z
```

This computes the reaction potential from charges (periodic), adds the linear ramp from the applied field to produce the total potential (`potential_total.xvg`), detects bulk water regions, and reports the applied voltage and slope-based estimate.

### Plain output without xmgrace headers

```bash
python potential.py -s topol.tpr -f traj.xtc -xvg none -center
```

## Choosing the number of slices

If `-sl` is not specified, the tool auto-estimates a value based on the box length and number of atoms (typically 50-1000 slices). Guidelines:

- **Too few slices** (e.g., 10-30): Poor spatial resolution, smooths out real features.
- **Too many slices** (e.g., >1000): Each slice contains very few atoms per frame, leading to noisy charge density.
- **Good range**: 200-500 slices for a typical ~10 nm box (slice width ~0.02-0.05 nm). More frames help compensate for finer slicing.

With the **Fourier method** (default), the result is stable across slice counts. With `-classical`, the result can vary dramatically with the number of slices (see below).

---

## The boundary artifact problem: why classical integration fails

This section describes a fundamental issue with the classical real-space double integration method used by `gmx potential` (and reproduced by `-classical` in this tool), and how the Fourier-space method solves it.

### The problem

When computing the electrostatic potential from an MD simulation, the standard recipe is:

1. Divide the simulation box into thin slices along one axis (typically Z)
2. Compute the charge density rho(z) in each slice
3. Integrate rho(z) twice to obtain the potential psi(z) via the Poisson equation

In the classical form (Eq. 2 in Gurtovenko 2009), the boundary conditions are psi(0) = 0 and E(0) = 0, and the integration proceeds from z = 0 to z = L (the box length):

```
E(z)   = (1/epsilon_0) * integral_0^z rho(z') dz'
psi(z) = -integral_0^z E(z') dz'
```

For a symmetric POPC bilayer (centered so the membrane is in the middle of the box and water is at the edges), this approach produces potential profiles that are **strongly dependent on the number of slices**, with large asymmetries that vary erratically:

```
Slices   Peak (V)   Asymmetry psi(L)-psi(0) (V)
  50      0.55       -0.08
 100      0.81       +0.47
 200      0.69       +0.17
 300      0.75       +0.25
 500      0.69       +0.13
 800      0.61       -0.00
1000      0.65       +0.06
```

For a system that is perfectly symmetric, psi(L) should equal psi(0). Instead, the potential drop across the box fluctuates between -0.08 and +0.47 V depending on how many slices are used. The peak potential varies between 0.55 and 0.81 V. This behavior is identical in `gmx potential`.

### Root cause: molecule splitting at the periodic boundary

The artifact originates from how atoms are assigned to slices when the system is periodic.

After centering (shifting so the bilayer center of mass is at z = L/2), each atom's z-coordinate is wrapped into [0, L) using modular arithmetic (`z % L`). This ensures all atoms fall within the box. However, **molecules that straddle the wrapping boundary at z = 0 / z = L get torn apart**: some atoms of the molecule end up near z = 0 and others near z = L.

For a typical POPC/water system, approximately 400 water molecules are split this way in each frame. Consider a TIP3P water molecule at the boundary:

- Before wrapping: O at z = -0.3 A, H1 at z = 0.5 A, H2 at z = 0.7 A (molecule is intact)
- After wrapping: O at z = L - 0.3, H1 at z = 0.5, H2 at z = 0.7 (molecule is torn apart)

The oxygen (charge -0.834e) is now at the right edge of the box, while the hydrogens (+0.417e each) remain at the left edge. Each split molecule creates an **artificial dipole spanning nearly the entire box**. While the total charge of each molecule is zero, its charges are now distributed across opposite ends of the integration path.

When the Poisson equation is integrated from z = 0 to z = L, this artificial dipole contributes a linear drift to the potential. The magnitude of this drift depends on:

1. **The number of split molecules** -- varies between frames as molecules diffuse
2. **The exact positions of split atoms within their slices** -- changes with slice width, hence with the number of slices

This is why the artifact varies erratically with the slice count: different slice widths cause the split charges to land in different bins, changing the effective artificial dipole moment and thus the potential drift.

### What doesn't fix it

Several intuitive approaches were tested and found to be ineffective:

**Making molecules whole before wrapping.** Even when using a preprocessed trajectory where molecules are intact (e.g., via `gmx trjconv -pbc whole`), the centering step followed by `z % L` wrapping re-splits molecules at the new boundary. The artifact is identical to using a raw trajectory.

**Wrapping by molecule center of mass.** Instead of wrapping individual atoms, one can shift entire molecules based on their center of mass. This keeps molecules intact but atoms of molecules near the boundary extend beyond [0, L). These out-of-range atoms must be handled somehow:

- **Clipping** (placing them in the edge slice): Concentrates charge at the boundary, making the artifact much worse. At 1000 slices, the potential reached 47 V.
- **Modulo index wrapping** (`index % n_slices`): Effectively re-splits the molecule -- the out-of-range atoms land on the opposite side. Results are identical to simple atom wrapping.

**Scaling coordinates to a fixed box size.** Gurtovenko & Vattulainen recommend scaling all coordinates to the initial box size to avoid effects from box size fluctuations under NPT pressure coupling. Testing showed this has negligible effect -- the artifact is dominated by molecule splitting, not box fluctuations. The asymmetry pattern was essentially unchanged with or without scaling.

**Placing the wrapping boundary through the membrane.** Moving the boundary from the water phase (where many molecules cross it) to the membrane center (where fewer cross). This made the artifact much worse, because lipid headgroups carry much larger partial charges than water, and splitting a lipid headgroup creates a correspondingly larger artificial dipole.

### Fix 1: Sachs et al. correction

Sachs et al. (J. Chem. Phys. 121, 10847, 2004) proposed an alternative form of the Poisson equation where the potential is constrained to be equal on both sides of the box: psi(0) = psi(L). This amounts to subtracting a linear correction term:

```
psi_Sachs(z) = psi_classical(z) - (z/L) * psi_classical(L)
```

Gurtovenko & Vattulainen (Eq. 4-11) showed that the difference between the two forms is proportional to the total dipole moment of the system. For electroneutral systems, the correction is:

```
Delta(z) = z / (epsilon_0 * V) * P_total
```

where P_total is the system dipole moment and V is the box volume.

The Sachs correction eliminates the linear drift and produces stable results:

```
Slices   Peak (V)   Asymmetry (V)
  50      0.56       -0.001
 100      0.58       +0.005
 200      0.61       +0.001
 300      0.63       +0.001
 500      0.62       +0.000
 800      0.61       -0.000
1000      0.62       +0.000
```

However, it has a fundamental limitation: **it removes ALL linear drift, including any real transmembrane potential difference**. For asymmetric bilayers (e.g., POPC/POPE), the transmembrane potential is a real physical quantity that should not be subtracted. The Sachs form is therefore only appropriate for symmetric systems.

### Fix 2: Dipole correction (subtracting the artifact only)

An alternative approach is to compute the artificial dipole moment caused by molecule splitting and subtract only that contribution:

1. For each frame, compute the dipole moment of the system with intact (unwrapped) molecules: P_whole = sum(q_i * z_i_whole)
2. Compute the dipole moment after wrapping: P_wrapped = sum(q_i * z_i_wrapped)
3. The artifact dipole is: P_artifact = P_wrapped - P_whole
4. Subtract the corresponding linear potential: Delta(z) = z / (epsilon_0 * V) * P_artifact

In principle, this should remove only the wrapping artifact while preserving the real transmembrane potential. In practice, testing showed this **overcorrects** -- the "whole" and "wrapped" dipole moments are not cleanly separable because the centering operation itself changes the dipole moment relative to the box origin. The corrected asymmetry was worse than uncorrected for many slice counts.

### Fix 3: Fourier-space integration (default)

The most robust solution is to solve the Poisson equation in Fourier (reciprocal) space, where periodicity is built in by construction.

The 1D Poisson equation in real space:

```
d^2 psi(z) / dz^2 = -rho(z) / epsilon_0
```

becomes, after Fourier transform:

```
-k^2 * psi_hat(k) = -rho_hat(k) / epsilon_0
```

which gives:

```
psi_hat(k) = rho_hat(k) / (k^2 * epsilon_0)    for k != 0
```

The k = 0 component corresponds to an arbitrary overall offset of the potential and is set to zero. The electric field is obtained similarly:

```
E_hat(k) = -i*k * psi_hat(k)
```

The inverse Fourier transform gives the real-space potential and field.

This approach has several key advantages:

1. **Inherently periodic.** The Fourier transform assumes the input signal is periodic. The potential automatically satisfies psi(0) = psi(L), with no need for corrections or special boundary conditions. There is no "starting point" for integration that could accumulate errors.

2. **Immune to the boundary artifact.** Split molecules create charge density features at both edges of the box. In the classical method, the double integration treats these as separated by the full box length, amplifying their effect. In Fourier space, periodicity means the edges are the same point -- the charge distribution from a split molecule is treated as if the molecule were whole, because the Fourier basis functions wrap around naturally.

3. **Insensitive to the number of slices.** Because there is no error accumulation from sequential integration, the result converges rapidly and varies minimally with slice count:

```
Slices   Peak (V)   Asymmetry (V)
  50      0.65       +0.003
 100      0.63       +0.003
 200      0.64       +0.001
 300      0.64       +0.001
 500      0.63       +0.000
 800      0.63       +0.000
1000      0.63       +0.000
```

The peak potential is consistent at ~0.63 V regardless of slice count, and the asymmetry is essentially zero (< 0.003 V). Compare this with the classical method's range of 0.55-0.81 V for the peak and up to 0.47 V of asymmetry.

### Summary of methods

| Method | Asymmetry | Peak stability | Works for asymmetric bilayers? |
|--------|-----------|---------------|-------------------------------|
| Classical (gmx potential) | up to 0.47 V | varies 0.55-0.81 V | yes, but with artifacts |
| Classical + Sachs | < 0.005 V | 0.56-0.63 V | no (removes real potential drop) |
| Classical + correct | < 0.004 V | 0.61-0.63 V | unclear (see discussion below) |
| Classical + dipole correction | overcorrects | overcorrects | no (unreliable) |
| **Fourier (default)** | **< 0.003 V** | **0.63-0.65 V** | **yes** |

The Fourier method is the default because it is the most robust, produces consistent results regardless of slice count, and works correctly for both symmetric and asymmetric systems.

---

## How it works

### Step 1: Charge density

The box is divided into `N` slices along the chosen axis. For each trajectory frame, every atom is assigned to a slice based on its coordinate. The partial charge (from the `.tpr` file) is added to that slice. The charge density is the total charge in each slice divided by the slice volume, averaged over all frames.

### Step 2: Solve the Poisson equation

**Fourier method (default):** The charge density is Fourier-transformed, divided by `k^2 * epsilon_0` in reciprocal space (skipping k=0), and inverse-transformed to obtain the potential. The electric field is computed analogously using `-i*k * psi_hat(k)`.

**Classical method (`-classical`):** The charge density is numerically integrated twice using the cumulative trapezoidal rule, with boundary conditions psi(0) = 0 and E(0) = 0. With `-sachs`, a linear correction term is subtracted to enforce psi(0) = psi(L).

Both methods use the `gmx potential` sign convention, where the dipole potential inside the membrane is positive relative to water.

## Known issues and considerations

- **Centering matters**: Without `-center`, bilayer center-of-mass fluctuations smear the charge density and produce asymmetric profiles. Always use `-center` for membrane systems.
- **Convergence**: Potential profiles require long trajectories (100+ ns for bilayers) to converge, particularly in the water region near lipid headgroups.
- **Dielectric constant**: The calculation uses vacuum permittivity (epsilon_r = 1). The output represents the bare electrostatic potential from the explicit charge distribution, not a dielectrically screened quantity.
- **Undulations**: For large membranes with significant undulations, the 1D slab decomposition may not accurately capture the local potential. Using small, planar bilayer patches avoids this issue.
- **Preprocessing**: It is recommended (but not strictly necessary) to preprocess the trajectory with `gmx trjconv -pbc whole` to make molecules whole. While the Fourier method handles split molecules much better than the classical method, intact molecules produce marginally cleaner charge density profiles.

## Output format

Output files use the GROMACS `.xvg` format, compatible with xmgrace/xmgr:

```
# This file was created by potential.py
# Electrostatic potential calculation using MDAnalysis
@    title "Electrostatic potential"
@    xaxis  label "Z (nm)"
@    yaxis  label "Potential (V)"
@TYPE xy
@ s0 legend "System"
  0.025000  0.000000
  0.075000  -0.012345
  ...
```

Use `-xvg none` for plain whitespace-separated columns without xmgrace headers.

## Differences from gmx potential

| Feature | `gmx potential` | `potential.py` |
|---------|----------------|----------------|
| Integration method | Classical double integration (artifacts) | Fourier (default) or classical (`-classical`) |
| Topology input | `.tpr` only | `.tpr` via MDAnalysis |
| Trajectory formats | `.xtc`, `.trr`, etc. | Any format MDAnalysis supports |
| Atom selection | GROMACS index groups | MDAnalysis selection syntax |
| `-correct` flag | Subtracts mean rho and mean E | Available via `-classical -correct` |
| Slice auto-estimation | No (default: 10) | Yes (50-1000 based on system size) |
| Sachs et al. form | Not available | Available via `-classical -sachs` |
| Fourier integration | Not available | Default method |
| Slice-count sensitivity | Severe | Negligible (Fourier) |
| Applied electric field | Not supported | `-efield` with water detection and voltage reporting |

The `-correct` flag from `gmx potential` is available via `-classical -correct`. See the detailed discussion below.

---

## The `-correct` flag: what it does and whether it's valid

### What gmx potential `-correct` actually does

Inspection of the GROMACS source code (`gmx_potential.cpp`) reveals that `-correct` performs **two mean subtractions** during the classical double integration:

1. **Before the first integration**: subtract the mean charge density from all slices that have nonzero charge density. This makes the effective total charge exactly zero, removing a constant drift term from the electric field.

2. **After the first integration**: subtract the mean electric field (again, only from slices with nonzero charge density). This removes a constant offset from E(z), which would otherwise produce a linear drift in the potential.

This is *not* the same as "zeroing net charge per slab" as sometimes described. It is a global correction applied to the entire charge density profile and then the entire electric field profile.

### Why it works

The key insight is that Step 2 (mean E subtraction) is doing essentially the same thing as the Sachs correction and the Fourier method's k=0 exclusion, but from a different angle.

In a periodic system, the electric field must also be periodic: E(0) = E(L). The classical integration starts with E(0) = 0, so it requires E(L) = 0 as well. If the total integrated charge density is not exactly zero (due to molecule splitting at the boundary, finite precision, or slight charge imbalance from the slicing), E(L) != 0, and this mismatch propagates as a linear drift in the potential.

Subtracting the mean electric field shifts E(z) so that its average is zero. For a periodic function, this is equivalent to enforcing that the integral of E over the full box vanishes — which is exactly the condition for the potential to be periodic (psi(0) = psi(L)).

Step 1 (mean rho subtraction) has a negligible effect in practice because the total system charge is already very close to zero. Its main effect is to eliminate tiny numerical imbalances from the binning.

### Test results

Testing on a symmetric POPC bilayer with the `-correct` flag produces results nearly identical to the Fourier method:

```
Classical + correct:
Slices   Peak (V)   Asymmetry (V)
  50      0.56       +0.002
 100      0.61       -0.003
 200      0.63       +0.001
 300      0.63       +0.001
 500      0.63       +0.000
 800      0.62       +0.000
1000      0.63       +0.000
```

Compare with Fourier (peak ~0.63 V, asymmetry < 0.003 V) and uncorrected classical (peak 0.55-0.81 V, asymmetry up to 0.47 V). The `-correct` flag effectively eliminates the slice-count sensitivity.

### Is it valid?

The `-correct` procedure is mathematically equivalent to enforcing periodic boundary conditions on the electric field, which is physically correct for a periodic system. In this sense, it is a valid fix.

However, there is an important subtlety for **asymmetric systems**:

- **Sachs correction** explicitly removes the linear potential drop across the box, destroying any real transmembrane potential difference.
- **`-correct`** subtracts the mean E field. For an asymmetric system with a real transmembrane potential, the mean E field is nonzero (it equals the potential drop divided by the box length). Subtracting it would remove the real transmembrane potential along with the artifact.

In practice, **`-correct` behaves identically to the Sachs correction for the linear component**: both remove the entire psi(L) - psi(0) drop. The difference is only in higher-order effects from the mean rho subtraction in Step 1 of `-correct`, which are negligible.

**Bottom line**: `-correct` is valid for symmetric systems and produces results equivalent to the Fourier method. For asymmetric bilayers with an intrinsic transmembrane potential arising from lipid asymmetry, `-correct` may partially remove that potential (like Sachs), since both subtract a global mean. However, for systems with an applied external electric field (see [below](#simulations-with-applied-electric-field)), `-correct` preserves the voltage-carrying slope in the reaction potential — the global mean E field that gets subtracted is dominated by the integration drift artifact, not the physical field structure. The Fourier method remains the recommended default because it handles all cases correctly and does not require any post-hoc corrections.

---

## Simulations with applied electric field

A common approach to simulate a transmembrane voltage in MD is to apply a constant external electric field E perpendicular to the membrane (GROMACS `electric-field-z` mdp parameter). The resulting voltage is V = E * L_z, where L_z is the box length in the field direction. This section explains what the tool computes in this scenario, how to interpret the results, and what the `-efield` flag does.

### Background: the constant electric field method

In GROMACS, a constant electric field E_z applies a force F = q_i * E_z to every charged particle. The applied voltage across the simulation box is:

```
V = E_z * L_z
```

This is equivalent to the influence of two infinite baths connected by an electromotive force (EMF), as shown by Roux (Biophys. J. 95, 4205, 2008). The constant field method has been validated extensively (Gumbart et al., BBA 1818, 294, 2012; Jensen et al., J. Gen. Physiol. 141, 619, 2013) and is widely used for studying ion channels, electroporation, and voltage-gated conformational changes.

### Decomposition of the total potential

The total electrostatic potential in a simulation with an applied field has two components (Jensen et al. Eq. 2):

```
Phi_total(z) = Phi_reaction(z) - E_applied * z
```

where:

- **Phi_reaction(z)** is the reaction potential arising from the explicit charges in the system (ions, water dipoles, lipid headgroups, protein residues). This is computed by the particle-mesh Ewald (PME) method during the simulation and is what our tool calculates from the charge density by solving the Poisson equation. **This potential is periodic**: Phi_reaction(0) = Phi_reaction(L).

- **-E_applied * z** is the linear ramp from the applied external field. This component has no source charges — it is imposed externally. **This breaks periodicity**: the total potential differs by V = E * L_z between the two ends of the box.

The charge density rho(z) that the tool computes from atom positions reflects the system's *response* to the applied field — the rearrangement of ions, water dipoles, and other charges. It does not contain the external field itself. Therefore, the Poisson equation solved by the tool yields only the reaction potential, which is periodic. The Fourier method remains the correct solver for this component.

### What happens in the bulk water

Bulk aqueous electrolyte behaves approximately as a conductor: it self-organizes to screen the applied field. In the water region:

```
E_total = E_reaction + E_applied ≈ 0
```

Therefore the reaction electric field in water is approximately:

```
E_reaction ≈ -E_applied
```

and the reaction potential develops a slope in the water region:

```
dPhi_reaction/dz ≈ E_applied    (in water)
```

The total potential, being the sum of the sloped reaction potential and the opposite-sloped applied ramp, is approximately flat in bulk water — as expected for a conductor. The entire voltage drop concentrates across the low-dielectric membrane region, which is the biologically meaningful quantity.

### Recovering the voltage from the charge density

A key question is whether the applied voltage V can be recovered purely from the charge density profile, without knowing E_applied a priori.

**From the electric field profile**: The average reaction electric field in the bulk water region should be approximately -E_applied. Reading this value and multiplying by L_z gives an estimate of V.

**From the potential slope**: The reaction potential has a slope ≈ E_applied in the water region. Fitting a line to the potential in the water phase and multiplying the slope by L_z gives V.

However, these are **approximate** estimates. Testing on a gramicidin A channel system (E_applied = 0.032 V/nm, expected V = 293 mV) showed:

```
Reaction potential slope in water (per region):
  Region 1 (left water):  slope = 0.0411 V/nm
  Region 2 (right water): slope = 0.0137 V/nm
  Average slope:           0.0274 V/nm
  -> V estimate:           251 mV (85.7% of expected 293 mV)

Average reaction E field in water: -0.0278 V/nm
Expected (-E_applied):             -0.0320 V/nm
Average total E field in water:     0.0042 V/nm (ideal: 0)
```

The ~85% recovery is consistent and reproducible. The underestimate reflects the fact that water is not a perfect conductor:

1. **Finite water phase**: The bulk water regions (~2 nm on each side) may not be thick enough for complete screening.
2. **Ionic current**: In a conducting channel (like gramicidin), steady-state ionic current requires a residual driving field in the water, so the screening is inherently incomplete.
3. **Finite system size**: Gumbart et al. (2012) emphasize that V = E * L_z is exact by construction, but extracting V from the potential profile introduces finite-size errors.

**The exact voltage is always V = E * L_z.** The slope-based estimate is useful as a consistency check but should not be used as the primary voltage determination.

### What the `-efield` flag does

When `-efield VALUE` is specified (where VALUE is E_z in V/nm, matching the GROMACS `electric-field-z` mdp parameter):

1. **Reports the exact voltage**: V = E_z * L_z, using the average box length from the trajectory.

2. **Writes the total potential**: Phi_total(z) = Phi_reaction(z) - E * z is written to the file specified by `-ot` (default: `potential_total.xvg`). This non-periodic profile shows the voltage drop across the membrane.

3. **Detects bulk water regions**: Water molecules are automatically identified (by residue name: SOL, TIP3, HOH, WAT, SPC) and their positions are binned during trajectory processing. Contiguous regions with high water density are identified as bulk water.

4. **Measures the reaction potential slope**: A linear fit is performed in each bulk water region independently (since different regions are at different absolute potential levels in the periodic reaction potential), and the average slope is reported.

5. **Reports voltage estimates**: Both the exact V = E * L_z and the slope-based estimate are printed, along with the recovery percentage. The average reaction electric field in water is also reported and compared to the expected -E_applied.

### Classical integration with `-correct` also works

Perhaps surprisingly, `-classical -correct` produces the same water slope and voltage estimate as the Fourier method. Testing on the gramicidin system shows identical slopes (0.0270 V/nm) for both methods.

The reason: `-correct` subtracts the global mean E field across all slices. This mean (-0.051 V/nm in the test system) reflects contributions from both the water regions (where E_reaction ≈ -0.029 V/nm) and the membrane region (where E_reaction ≈ +0.017 V/nm, reflecting the large headgroup dipole fields). Because these contributions partially cancel, the global mean is dominated by the integration drift artifact, not by the physical field from the applied voltage. Subtracting it removes only the artifact, preserving the spatial structure — exactly as the Fourier method does by setting k=0 to zero.

```
                           Uncorrected    After -correct    Fourier
Mean E, all slices:        -0.0512         0.0000           0.0000  V/nm
Mean E, left water:        -0.0826        -0.0314          -0.0314  V/nm
Mean E, right water:       -0.0775        -0.0263          -0.0263  V/nm
Mean E, membrane:          -0.0338        +0.0174          +0.0174  V/nm
```

Both corrected methods give identical E field profiles. The `-sachs` correction, by contrast, enforces psi(0) = psi(L) by subtracting a linear term, which would remove the physical slope and should not be used with applied fields.

### Practical considerations

- **Always use `-center`** when analyzing systems with applied fields. Without centering, membrane drift smears the charge density and corrupts the slope measurement.

- **The total potential (`-ot` output) is not periodic.** It shows a linear ramp across the box with the voltage drop concentrated at the membrane. This is the physically meaningful potential profile.

- **The reaction potential (`-o` output) is periodic.** It shows the intrinsic membrane dipole potential plus the system's response to the applied field. The slope in the water region reflects the screening of the applied field.

- **Both Fourier and `-classical -correct` work correctly** with applied fields. Both methods remove only the unphysical DC offset (k=0 mode) from the electric field while preserving all spatial structure, including the slope in water that carries the voltage information. The `-correct` procedure subtracts the global mean E field across all slices — this mean (-0.051 V/nm in the gramicidin test) is NOT the same as the E field in water (-0.029 V/nm), because the membrane region has a very different E field that pulls the average away. Therefore, the physical slope is preserved. The `-sachs` correction should still be avoided with applied fields, as it explicitly enforces psi(0) = psi(L), which would remove the linear component of the reaction potential that encodes the voltage.

### Why the tool computes the reaction potential, not the total potential

A potentially confusing aspect of electrostatic potential calculations from MD simulations — both in `gmx potential` and in this tool — is that what gets computed is **the reaction potential**, not the total electrostatic potential. This distinction is invisible in simulations without an applied field, which is why it has historically caused confusion.

**What the tool does:** It bins atomic charges into slices, constructs the charge density rho(z), and solves the Poisson equation:

```
d²Phi/dz² = -rho(z) / epsilon_0
```

The charge density rho(z) contains only the explicit charges in the system: ions, water partial charges, lipid headgroup charges, protein residues, etc. The resulting potential Phi(z) is the **reaction potential** — the electrostatic potential generated by the system's own charges in response to whatever conditions are imposed.

**Why not the total potential?** When an external electric field is applied in GROMACS (`electric-field-z`), the field acts as a force F = q * E_z on every particle, but it has **no source charges**. It is an external boundary condition, not a charge distribution. Since the applied field contributes no charges to rho(z), it does not appear in the Poisson equation, and the tool cannot "see" it. The total potential is:

```
Phi_total(z) = Phi_reaction(z) - E_applied * z
```

The linear ramp from the applied field must be added manually (which is what the `-efield` flag does).

**Why this hasn't caused problems historically:** In the vast majority of membrane simulations — symmetric bilayers, asymmetric bilayers, systems with ion gradients — there is no applied external field. In these cases, Phi_reaction(z) *is* the total potential. The distinction between reaction and total potential simply doesn't arise, and the output of `gmx potential` is the complete answer.

**Where confusion can arise:** When researchers apply an external field (e.g., to study voltage-gated channels or electroporation) and then run `gmx potential` on the trajectory, they get the reaction potential — which is periodic and shows no net voltage drop across the box. This can be mistaken for a calculation error or a sign that the applied voltage had no effect. In fact, the tool is working correctly; the voltage drop is carried entirely by the applied field ramp, which is not part of the charge density. The reaction potential shows the system's *response* to the applied voltage: screening in the water phase (visible as a slope) and the intrinsic membrane dipole potential.

**To obtain the total potential** in simulations with an applied field, use `-efield VALUE` where VALUE matches the `electric-field-z` parameter from the GROMACS mdp file. The tool will then write the total potential (reaction + applied ramp) to `potential_total.xvg`, which shows the non-periodic voltage profile with the full voltage drop across the membrane.

### Why is the reaction potential periodic?

At first glance, it seems paradoxical that the reaction potential Phi_reaction(z) is periodic even when a directional electric field is applied. One might expect charges to accumulate on one side of the membrane — as in a real capacitor — breaking the symmetry and making the charge density (and hence the potential) non-periodic. But this does not happen, for a fundamental reason: **periodic boundary conditions prevent true charge accumulation**.

**There is only one water phase.** In a membrane simulation with PBC, the water region above the membrane and the water region below it are not separate compartments — they are the same continuous water phase, connected through the periodic boundary. A water molecule (or ion) that exits the box at z = L_z immediately re-enters at z = 0. There are no edges, walls, or capacitor plates where charge can build up.

**Charge density inherits the periodicity of the box.** Because every atom's position is defined modulo L_z, the ensemble-averaged charge density satisfies rho(z + L_z) = rho(z) by construction. No matter how strong the applied field, the charge density profile repeats with period L_z. The system responds to the field through *local* polarization — water dipoles align preferentially, ions redistribute slightly within the periodic cell — but these rearrangements are periodic.

**Poisson's equation preserves periodicity.** The reaction potential is obtained by solving the Poisson equation with the periodic charge density as input. If rho(z) is periodic, Phi_reaction(z) is also periodic (the Fourier method enforces this explicitly; the classical method with `-correct` achieves it by removing the DC drift). The reaction potential develops a slope in the water regions (reflecting the screening of the applied field) but this slope is compensated by the opposite response in the membrane region, so that the potential returns to its starting value over one full period.

**How does the system sustain a voltage without charge accumulation?** This is the key insight from Roux (2008). In PBC, the applied voltage cannot be represented by accumulated surface charge (there are no surfaces). Instead, it enters as an external parameter — the EMF in Roux's formulation — that adds a force F = q * E_z to every charged particle. The system's response is captured by the *displacement charge* Q_d = sum(q_i * z_i) / L_z, which measures how much charge has been displaced along z relative to the periodic cell, without ever breaking periodicity. The displacement charge is analogous to the charge on a capacitor plate, but it arises from the collective shift of all charges within the periodic cell, not from accumulation at a boundary.

**In summary:** PBC means the simulation box has no boundaries — only a periodic unit cell. Charges cannot accumulate at edges that don't exist. The charge density and reaction potential are therefore necessarily periodic. The entire voltage drop is encoded in the external field ramp (-E * z), which must be added manually to obtain the non-periodic total potential.

### Special case: pure water under applied electric field

Testing the tool on a pure water box with and without an applied electric field reveals a counterintuitive but physically correct result: **the reaction potential is identical in both cases**, and the total potential under the applied field is a perfectly linear ramp.

**Why the charge density doesn't change:** The applied field does orient water dipoles — each molecule tilts slightly so its dipole aligns with E_z. However, this polarization is *spatially uniform*: every slice of the water box has the same degree of alignment. Uniform polarization produces no charge density variation:

```
rho_bound = -dP/dz = 0    (when P is constant)
```

Microscopically: if all hydrogens shift slightly "up", each slice gains H charge from molecules below but loses the same amount to the slice above. For a uniform density of molecular centers, these contributions cancel exactly in every slice. The net dipole moment is real, but it is invisible in rho(z).

**Why there is no screening:** In a finite dielectric slab (or at a membrane interface), the edges of the polarized region carry bound surface charges (sigma_b = P dot n) that create a depolarization field screening E_applied by a factor of 1/epsilon_r. In PBC, there are no surfaces — the "+sigma" at z = L_z wraps around and cancels the "-sigma" at z = 0. Without bound surface charges, there is no screening:

| System | E_reaction in water | E_total in water | Screening? |
|--------|-------------------|------------------|------------|
| Pure water (PBC) | 0 | E_applied | No |
| Water + membrane | ~ -E_applied | ~ 0 | Yes |

This means that in pure water with PBC, the total electric field equals the full applied field everywhere. The dielectric constant of water is irrelevant — epsilon_r only matters when interfaces are present for bound charges to accumulate on. The voltage V = E * L_z is exact regardless of epsilon_r, which is a fundamental property of the constant-field method in PBC (Roux, 2008).

**What this means in practice:** The perfectly linear total potential ramp in pure water is NOT because water is "a good conductor that screens the field" (leaving only the ramp). It is the opposite — water *cannot* screen the field at all in homogeneous PBC, so the reaction potential is flat and the total potential is the bare applied ramp. Screening requires interfaces (like a membrane) where the polarization changes spatially, creating bound charges that generate the reaction field.

**Connection to uncompensated dipole artifacts:** This physics also explains the PBC-induced water ordering artifact described by Kopec & Gapsys (2026): when a membrane has an uncompensated electric dipole (from asymmetric charged lipids or an embedded protein), PBC constrains the voltage across the box to zero. The system compensates by creating a counter-field — water dipoles order and ions redistribute. This ordering is an artifact of the periodic boundary, not a physical response. Adding salt (mobile ions) alleviates it by providing charge carriers that can redistribute to screen the membrane dipole without forcing bulk water to order.

### References

- B. Roux, "The membrane potential and its representation by a constant electric field in computer simulations," *Biophys. J.* **95**, 4205-4216 (2008).
- J. Gumbart, F. Khalili-Araghi, M. Sotomayor, B. Roux, "Constant electric field simulations of the membrane potential illustrated with simple systems," *BBA - Biomembranes* **1818**, 294-302 (2012).
- M.O. Jensen, V. Jogini, M.P. Eastwood, D.E. Shaw, "Atomic-level simulation of current-voltage relationships in single-file ion channels," *J. Gen. Physiol.* **141**, 619-632 (2013).
- W. Kopec, V. Gapsys, "Periodic boundaries in molecular dynamics simulations: why do we need salt?," (2026).
