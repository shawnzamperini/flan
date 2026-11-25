# HERMES Backgrounds

Flan can follow impurities in a HERMES plasma background, set with
`bkg_source = "hermes"`. HERMES here means a Hermes-3 turbulence run on
a hypnotoad single-null grid, either from JAX-Hermes3 (`hermes_format =
"jax"`, the default) or from BOUT++ Hermes-3 (`hermes_format = "bout"`).
The two write their output differently, but Flan builds the same
background from either. This page explains
how HERMES backgrounds work in Flan. An example input is
`regression/hermes_example.cpp`, and the options are listed under
[HERMES](../input-options.md#hermes) on the Input Options page. How the
coupling was built, and the assumptions it makes about HERMES, are on the
[HERMES Coupling Plan](../developer/hermes-coupling.md) page.

!!! note "GPUs"
    HERMES backgrounds run on GPUs (`use_gpu = "cuda"`) as well as the CPU,
    field-aligned interpolation, toroidal wedge, thermal force and target
    sputtering included, except for self-sputtering, which stops, and the
    sheath field, which a GPU run goes without. These are CPU only so
    far.

## Files

Both formats need the run's hypnotoad grid (`hermes_grid_file`,
netCDF): the cell corners, the magnetic field (`Brxy`, `Bzxy`, `Btxy`),
the field-line shift `zShift` and `ShiftAngle`, and the separatrix and
X-point indices (`ixseps1`, `jyseps1_1`, `jyseps2_2`).

In both, the output times to load are picked by index, in time order,
with `hermes_snap_start`, `hermes_snap_end` and `hermes_snap_stride`;
negative indices count from the end. At least two must be loaded.

### JAX-Hermes3

The snapshots (`snap*.npz` in `hermes_snap_dir`), one per output time,
each with `Ne`, `Pe`, `Pi`, `NVi`, `phi` and `t` on the (x, y, z) cell
centers, in name (time) order.

The snapshots are normalized, and the normalization is not in them. The
`hermes_*` options must be set from the case's input JSON:

| JSON (`pars`) | Option |
|---|---|
| `norm` (Nnorm, Tnorm, Bnorm) | `hermes_nnorm`, `hermes_tnorm`, `hermes_bnorm` |
| `average_atomic_mass` (or the ion's `AA`) | `hermes_ion_mass_amu` |
| `zperiod` | `hermes_zperiod` |
| `kappa_coefficient_e`, `kappa_coefficient_i` | `hermes_kappa_coef_e`, `hermes_kappa_coef_i` |

### BOUT++ Hermes-3

`hermes_snap_dir` is the run directory, with the dump files
`BOUT.dmp.0.nc`, `BOUT.dmp.1.nc`, ..., one per MPI process. Each holds
its process's block of cells, with guard cells, at every output time.
Flan reads the layout (`NXPE`, `NYPE`, `MXSUB`, `MYSUB`, `MXG`, `MYG`,
`MZ`), the output times (`t_array`), the normalization (`Nnorm`,
`Tnorm`, `Bnorm`) and `zperiod` from the dumps, and assembles the blocks
into the global grid by each file's `PE_XIND` and `PE_YIND`. A run split
in z (`NZPE` > 1) isn't supported. One run directory makes one
background: restarts written to another directory can't be joined.

The variables read are `Ne`, `Pe`, `P<ion>`, `phi`, the ion flow, and
for an electromagnetic run `Apar`. `<ion>` is the name of the main ion
species in the run's input, set with `hermes_ion_species` (`"i"` by
default, so `Pi`, `Vi`). The ion flow is `V<ion>` if it was written,
otherwise it comes from `NV<ion>`: in an electromagnetic run that is the
canonical momentum AA N V + N Apar, so the flow is (NV - N Apar) /
(AA N).

The dumps don't have the ion mass or the conduction coefficients. Set
them from the run's `BOUT.settings` (which
`python/check_hermes_conventions.py --format bout --run-dir` prints):

| `BOUT.settings` | Option |
|---|---|
| `AA` of the main ion species | `hermes_ion_mass_amu` |
| `kappa_coefficient` of the electrons and ions | `hermes_kappa_coef_e`, `hermes_kappa_coef_i` |
| `kappa_limit_alpha` (if set) | `hermes_kappa_limit_alpha` |

`hermes_nnorm`, `hermes_tnorm`, `hermes_bnorm` and `hermes_zperiod` are
not used.

**Inductive E∥.** In an electromagnetic run E∥ = -∇∥φ - ∂A∥/∂t. With
`hermes_inductive_epar = "on"` (the default), Flan adds -∂A∥/∂t when the
dumps have `Apar`. For each loaded output time it reads `Apar` at the
output times just before and after, and keeps the centered difference
(one sided at the run's first and last output times), so it doesn't
depend on `hermes_snap_stride`. It is only as good as the output
interval is short compared with how fast `Apar` changes. Electrostatic
runs, and JAX-Hermes3, have no `Apar` and use E∥ = -∇∥φ.

**Reading time and memory.** The dump files hold every output time, so
each file is opened once and all the loaded times are read from it
together. With many files the reading takes minutes, mostly from the
number of files; Flan prints its progress. The background holds five 3D
fields per loaded time (more once E, the flows and the heat fluxes are
added), and Flan prints the size of the plasma before it starts reading.
Use `hermes_snap_stride` to cover a long time window in less memory.

### The plasma

The plasma is n_e = Ne Nnorm, T_e = (Pe / Ne) Tnorm and T_i = (Pi / Ne)
Tnorm, with one ion species of charge 1 (n_i = n_e). The parallel ion
flow (NVi / (AA Ne) for JAX-Hermes3, see above for BOUT++) is in units of
sqrt(e Tnorm / m_p), the potential is phi
Tnorm and the time is t m_p / (e Bnorm). Cells with Ne <= 0 get zero
temperatures and flow. The electric field is minus the gradient of the
potential, and the background ion flow is the parallel flow plus the
E x B and grad-B drifts, as for the other backgrounds.

## Zones and coordinates

The grid is split into six zones at the separatrix and the X-point:

| Zone | Region |
|---|---|
| 1 | Inner private flux region (PFR) |
| 2 | Core |
| 3 | Outer PFR |
| 4 | Inner divertor leg (SOL) |
| 5 | Main SOL |
| 6 | Outer divertor leg (SOL) |

In every zone, x and y are the grid's global cell indices and z is the
toroidal angle in radians:

- x counts the grid's two guard cells, so the separatrix is at
  x = `ixseps1`, and x = 20.5 is the middle of cell 20.
- y runs from the inner target (y = 0) to the outer one, in the grid's
  order.
- z covers 1 / `hermes_zperiod` of the torus, from -dz/2 to
  2π / `hermes_zperiod` - dz/2 with dz the toroidal cell width, as in
  BOUT++.

The starting position of particles (`imp_xstart_opt`, etc.) is given in
these coordinates, in the zone set by `imp_zone_start_val`. (R, z, Z) is
taken as right handed.

## Boundaries

- The y ends of the legs and PFRs are the divertor targets, which are
  wall faces.
- The outer edge of the SOL and the PFRs (the grid's last x cells) is an
  absorbing face.
- The core's inner edge is a core face: a particle that crosses it is put
  back at a random place on it.
- In z the domain is periodic. When it is a wedge of the torus
  (`hermes_zperiod` > 1), a particle that leaves through one end comes
  back in at the other, a wedge's width around the torus, and its
  velocity is rotated with it.

## How the background is built

**Along the field line.** On a hypnotoad grid the field line through a
cell moves many toroidal cells between neighboring y cells, so in a
turbulent run the plasma at the same z in the next y cell is unrelated.
Flan therefore interpolates the plasma in y along the field line, using
`zShift`, and takes y derivatives along it, shifting the neighboring
planes in z with a Fourier transform as BOUT++ does. The parallel
gradient, used for E∥ and the heat fluxes, is HERMES's own:
∇∥f = (B·∇y / B) ∂f/∂y along the field line. `hermes_field_aligned =
"off"` does all of this at fixed z instead, and is only meant for
comparison.

**Parallel heat fluxes.** The runs' output doesn't have them, so for the
thermal force Flan calculates them the way Hermes-3 does
(`hermes_heat_flux = "conduction"`):
q∥ = -κ∥ ∇∥T for electrons and ions, with κ∥ = c n e² T τ / m (T in eV),
τ = 1 / Σν from all of the species' Braginskii collisions (e-e and e-i
for electrons, i-i and i-e for ions) and c = `hermes_kappa_coef_e` or
`_i`. `hermes_kappa_limit_alpha` turns on Hermes-3's flux limiter. The
heat fluxes are positive along B and saved with the background as
`q_par_e` and `q_par_i` (W/m²). As in Hermes-3, nothing is conducted
through a target face (the sheath carries that energy), so the q∥ of a
cell at a target is half the flux through its face with the next cell.
Without a flux limiter, the ion heat flux in a hot, weakly collisional
SOL can exceed n T v_th, beyond what the heat-flux-skewed sampling of the
thermal force can represent; there the skew is capped (|ε| ≤ 1), and the
run's output reports how often.

**Ghost cells** two cells deep are filled from the neighboring zones,
as for SOLEDGE3X. At the X-point, a corner ghost cell takes the average
of the cells around it, which lie on different field lines, so
interpolation in the cells touching the X-point is approximate.

## Sputtered source

With `imp_source_opt = "sputter"` the particles are the atoms the
background ions sputter off the divertor targets (see
[Sputtering](../physics/sputtering.md)). The grid ends at the targets, so
the cell faces on them are the sputtering surface, each with its own
normal and area. There's no wall mask or wall outline to correct them
with, and `sput_max_wall_dist` isn't used. The ion flux onto a target face
is the parallel flow through it, n u∥ b·n, plus the E×B drift into it,
from the plasma in the cell next to it. An atom that comes back to a target
is put on the face element it crossed, so the output has the gross and net
erosion along each target. Only the targets sputter: the outer edges of
the SOL and the PFRs are absorbing faces, not walls.

The domain is the run's toroidal wedge, so the sputtered rate, and the
weight of each particle, is for that wedge: 1 / `hermes_zperiod` of the
whole torus's. Densities don't change, since the wedge's volume shrinks
in proportion.

## Checking a background

`python/check_hermes_conventions.py` checks a run's files before they
are used: the cell corners, the handedness, the field-line shift, the
flow into the targets and the particle and energy balance.
`python/check_hermes_bkg.py` checks the background as Flan built it,
from a Flan output file: the parallel ion flow into each target against
the sound speed, the direction of the ion grad-B drift relative to the
X-point (`--expect toward` for the forward field, `away` for the
reversed), and, with `--content`, where the impurities are.

## Output

The output file has the same multi-zone layout as for SOLEDGE3X
backgrounds (see [Output](soledge3x.md#output)), without the wall mask
`chi`. `python/plot_multizone.py` plots the impurity density, or any
background quantity with `--var` (e.g., `Te` or `q_par_e`), in the
(R, Z) plane, for one toroidal plane of the wedge or averaged over them.

Each loaded output time holds the full 3D background, so memory grows
with the number of times loaded. `hermes_snap_stride` loads every n-th
one, to cover a long time window in less memory. The times are a stretch of a
turbulent run, so `tbound_type = "periodic"` loops over them and treats
the turbulence as statistically steady.

## Limits

- One ion species of charge 1 and quasineutrality (n_i = n_e). A run with
  several ion species would need their densities and collisions.
- No neutrals: the impurities see only the ions and electrons.
- The background is used output time by output time, interpolated linearly in
  time, with no time averaging.
- BOUT++ runs: one run directory per background, and no split in z
  (`NZPE` = 1).
- On GPUs, no sputtering or sheath yet (see the note at the top).
