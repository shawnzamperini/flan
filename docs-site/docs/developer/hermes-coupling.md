# HERMES Coupling Plan

This page is the development plan for reading HERMES plasma backgrounds
into Flan, the same way SOLEDGE3X backgrounds are read
(`bkg_source = "hermes"`). The first version reads output from
JAX-Hermes3, a Python/JAX reimplementation of Hermes-3. The plan for a
second reader, for BOUT++ Hermes-3 dump files, is in
[BOUT++ Hermes-3 reader](#bout-hermes-3-reader). The last section is a
short summary for the HERMES-3 developers of how the coupling works and
what it assumes.

## The data

A JAX-Hermes3 run has two parts:

- **Grid**: a hypnotoad netCDF file with `Rxy`, `Zxy`, the four cell
  corners, `Brxy`, `Bzxy`, `Btxy`, the metric, `zShift`, `ShiftAngle` and
  the separatrix indices. It has 2 x guard cells on each side and no y
  guard cells.
- **Snapshots**: one uncompressed `.npz` file per time holding `Ne`, `Pe`,
  `Pi`, `NVi`, `NVe`, `phi` and `Vort` on the (x, y, z) cell centers,
  normalized by Nnorm, Tnorm and Bnorm, plus the time `t`. Each case also
  has a JSON file with the normalization, the species masses, the
  conduction coefficients and `zperiod`. The domain is 1/`zperiod` of the
  torus.

A BOUT++ Hermes-3 run instead writes one netCDF file per MPI process
(`BOUT.dmp.N.nc`) holding every output time, with the normalization in
each file. The variables are named after the species in the input
(`Ne`, `Pe` for electrons `e`; `Ni`, `Pi`, `NVi`, `Vi` for ions `i`, or
`Nd+`, `Pd+`, ... for `d+`). The grid file is the same, but the plasma
needs more care than a change of file format (see
[BOUT++ Hermes-3 reader](#bout-hermes-3-reader)).

## Zones

A single-null BOUT++ grid splits into 6 Flan zones that share whole faces:
the inner private flux region (PFR), the core and the outer PFR inside the
separatrix (x < `ixseps1`), and the inner leg, main SOL and outer leg
outside it. The y ranges are split at `jyseps1_1` and `jyseps2_2`. The y
ends of the legs are the targets (`face_wall`), the inner x boundary of
the core is `face_core`, and the inner x boundary of the PFR and the outer
x boundary of the SOL are absorbing. The core's last y cell connects to
its first, and the inner leg's last y cell to the outer leg's first, in
the PFR.

## Conventions

These were checked on the data with `python/check_hermes_conventions.py`,
which should be run on any new grid or case:

- **Corners**: `Rxy_corners` and `Zxy_corners` are the (-x, -y) corner of
  each cell; in the other corner names "right" is +x and "upper" is +y.
  Cells on either side of the branch cuts share their corners.
- **Handedness**: (R, z, Z) is right handed, with z the toroidal angle
  increasing counterclockwise seen from above. The grid's Jacobian `J`
  and curvature `bxcvx` both agree with this. So a point is at
  (X, Y, Z) = (R cos z, R sin z, Z), unlike SOLEDGE3X, where Y has the
  opposite sign.
- **Field direction**: the poloidal field (`Brxy`, `Bzxy`) points along
  +y in every cell (hypnotoad orients y that way), so +y runs along the
  poloidal part of b. The sign of `Btxy` sets the toroidal direction of
  the field.
- **Parallel flow**: V = NV / (AA N), positive toward +y. It points into
  both targets, as the sheath boundary requires.
- **Toroidal coordinate**: z is the lab-frame toroidal angle (shifted
  metric). Along a field line, going from plane y to y + 1, z changes by
  `zShift(y + 1) - zShift(y)`, plus `ShiftAngle` across the core branch
  cut. Nothing has to be twisted where zones meet.
- **Grid `chi`**: a poloidal coordinate, not a wall mask. It isn't used.

## Phases

**Phase 0: conventions.** Done. `python/check_hermes_conventions.py`
checks the corners, the handedness and the grad-B drift direction, the
field alignment and the sign of the flow, and compares the particle and
energy fluxes reaching the targets with the sources. Open: confirm with
the JAX-Hermes3 developers that its conduction is the same as Hermes-3's
(below).

**Phase 1: refactor, no change in behavior.** Done. The code that isn't
specific to SOLEDGE3X (geometry, B and grad B, E and the flows, and the
derivatives they use) is in `axisym_mesh.cpp`, working from an
`Axisym::Mesh` in SI units, and the SOLEDGE3X reader converts its mesh to
one (`S3X::to_axisym`). The toroidal handedness is a property of the
background (`Background::get_phi_sign`, broadcast with it) instead of the
constant `S3X::phi_sign`, and sputtering uses it. A SOLEDGE3X background
read before and after the change was the same bit for bit, and
`test_axisym_mesh_cpu` tests the shared code on a mesh where everything
is known analytically, for both handednesses.

**Phase 2: readers and options.** Done.

- `npz.h/.cpp` reads uncompressed `.npz` files, taking the entry sizes
  from the zip's central directory, since numpy writes zip64 local
  headers without them. Compressed or Fortran-order files stop the run.
  Tested in `test_npz_cpu`.
- `read_hermes.cpp` reads the hypnotoad grid into the 6 zones (each with
  the grid's global cell indices as x and y), builds the
  `Axisym::Mesh` from the cell corners and B, and reads the snapshots:
  n_e = Ne Nnorm, T_e = (Pe / Ne) Tnorm, T_i = (Pi / Ne) Tnorm,
  phi Tnorm, and the flow along B, NVi / (AA Ne) c_s0, with its sign
  from the direction of the poloidal field relative to +y. Tested in
  `test_read_hermes_cpu` on a small grid and snapshots the test writes.
- Options `bkg_source = "hermes"` and `hermes_*` (see
  [Input Options](../input-options.md)). `hermes_format = "bout"` is
  accepted but stops with "not implemented yet".
- Known difference: Flan puts each cell's center at the average of its
  four corners, as for SOLEDGE3X, while hypnotoad's `Rxy`, `Zxy` (where
  B is given) lie on the curved flux surface. They differ by a few
  percent of a cell typically, up to about a fifth of a cell in the
  largest outer SOL cells.

**Phase 3: field-aligned interpolation and gradients.** Done. On a
shifted-metric grid a field line moves many toroidal cells between
neighboring y planes (typically tens of cells, and several turns of a
toroidal wedge near the X-point), so interpolating or
differentiating in y at fixed z mixes unrelated turbulent structures.

- `Zones` holds an optional `zShift(x, y)` table per zone, including the
  ghost cells, with `z_period` and `wrap_z`. The HERMES reader fills it
  from the grid: a ghost cell takes the zShift of the cell its plasma
  comes from, plus `ShiftAngle` where the field line crosses the core's
  branch cut. Zones without a table (slab, Gkeyll, SOLEDGE3X) are
  unchanged.
- Interpolation (`Utilities::interp_aligned`): in each of the two x
  columns, the field line through the point is followed to the two
  bracketing y planes, where the field is interpolated linearly in z
  (wrapped into the wedge). x stays at fixed z, as the shifted metric
  intends. n, T_e, T_i and q_par use it; B and the geometry don't, since
  they don't vary along z the way the turbulence does.
- Vectors (u, and E in the Boris push) are interpolated the same way,
  with each sample's (X, Y) part rotated from its cell's toroidal angle
  to the particle's. A sample taken a radian away would otherwise point a
  radian off. The periodic ghost cells at the ends of the wedge hold
  their vectors already rotated to their own angle (Phase 4).
- y derivatives (`Axisym::aligned_dfdy`): the neighboring planes are
  shifted onto the field line spectrally (`shift_periodic`, as BOUT++'s
  shifted metric does), then differenced, giving df/dy along the field
  line. The gradient uses df/dy at fixed z = df/dy|_ζ - (dzShift/dy) df/dz.
  The parallel gradient is taken from df/dy|_ζ alone,
  `grad_par` = (B·∇y / B) df/dy|_ζ, which is HERMES's own ∇∥. The
  parallel part of E is replaced by it. Flan's geometry (from the cell
  corners) and `zShift` agree to about 1% per cell, which is plenty for
  E⊥. But b·∇ of the full gradient then picks up about 1% of the much
  larger gradient across the field, which was as large as E∥ itself.
- Checks: unit tests in `test_field_aligned_cpu` on the test grid with a
  plasma constant along field lines (exact y derivative, no E∥, vectors
  rotated exactly, zShift agreeing with B and the basis). On HERMES
  output, E∥ is checked against a separate calculation of -∇∥φ along
  the field lines and against the electron parallel force balance,
  -∇∥p_e / (e n_e), which a nearly adiabatic edge plasma should follow
  closely. The fixed-z gradient fails both badly. Following the field
  line makes the particle loop roughly 20% slower on HERMES backgrounds.
- Corner ghost cells at the X-point have several sources on different
  field lines and get their average, as before, so interpolation within a
  cell of the X-point is approximate.

**Phase 4: toroidal wedge.** Done. On SOLEDGE3X and HERMES grids z is
the toroidal angle, and the domain may be a wedge of the torus (2π /
`zperiod` for HERMES, the span of `phicorners` for SOLEDGE3X). Its two z
faces are then the same surface a wedge's width apart, while Flan's
vectors are Cartesian, so they have to turn with the toroidal angle.

- `Zones::set_toroidal` marks a zone's z as the toroidal angle (both
  readers set it for every zone). `is_wedge` is true for a toroidal zone
  that's periodic in z with a period other than 2π.
- Ghost cells: `Zones::fill_ghosts(vX, vY, phi_sign)` rotates each
  source's (X, Y) part about Z by phi_sign (z_ghost - z_source) before
  copying, so a ghost past a z face holds the vector at its own angle.
  It's used for the cell centers, the reciprocal basis vectors, B,
  grad B, E and u. Without a wedge the rotation would be by 0 or 2π, so
  the components are copied as before, bit for bit.
- The field-aligned vector interpolation of Phase 3 rotates each sample
  from its cell's own angle (past the zone's range for a ghost), now
  that the ghosts hold rotated vectors.
- Particles: when a particle wraps in z on a wedge,
  `Boundary::check_bounds_cpu` rotates (vX, vY) by phi_sign times the
  change in z. A particle put back in at the core keeps its velocity, as
  before.
- The GPU doesn't rotate yet, so a wedge background stops with an error
  there (only single-zone backgrounds run on the GPU so far anyway).
- Checks: `test_wedge_cpu`, on the SOLEDGE3X test case as a quarter
  torus and as the whole torus: the ghost vectors match the
  axisymmetric field at their own angle, B interpolated on the two z
  faces differs by exactly the wedge's rotation, a wrapping particle's
  velocity is turned by the wedge's width, and the whole torus is
  unchanged.

**Phase 5: thermal force.** Done. The snapshots have no heat fluxes, so
Flan calculates them the way the run did, for electrons and ions:
q_par = -kappa_par grad_par T, positive along B like the SOLEDGE3X ones.

- `Braginskii` (`braginskii.cpp`): the collision frequencies and Coulomb
  logarithms of Hermes-3's `braginskii_collisions` and the conductivity
  of `braginskii_conduction`, kappa = c n e^2 T tau / m (T in eV), in SI
  units. tau = 1 / Σν uses Hermes-3's default "multispecies" set: e-e and
  e-i for electrons, i-i and i-e for ions. The optional flux limiter is
  included but off by default, as in Hermes-3.
- `Axisym::calc_conductive_heat_fluxes` takes grad_par from the
  field-aligned derivative of Phase 3 (not b·∇T, for the reason given
  there). The plasma is one ion species with n_i = n_e / Z; a run with
  several ion species would need their densities and collisions. The y
  derivative is one sided in the last cell before a target, so q_par
  there is approximate.
- Options `hermes_heat_flux`, `hermes_kappa_coef_e`, `hermes_kappa_coef_i`
  and `hermes_kappa_limit_alpha`. The heat fluxes are saved with the
  background and used by the thermal force as for SOLEDGE3X.
- Checks: `test_braginskii_cpu` (the frequencies against the NRL
  Formulary's Braginskii collision times, Spitzer's kappa_e, the Coulomb
  logarithm's ranges and floors, the limiter) and tests in
  `test_field_aligned_cpu` (the heat flux for T rising along the field
  line, none for T that only varies across it, the option). On HERMES
  output, Flan's q_par was checked against a separate calculation with
  BOUT++'s Grad_par from the grid's metric, and the conducted and
  convected power into the legs against the input power.

**Phase 6: validation and output.** Done.

- `regression/hermes_example.cpp`: tungsten from a point source just
  outside the separatrix at the outer midplane, with the thermal force,
  ionization and recombination, and each `hermes_*` option's source in
  the case's JSON.
- `hermes_field_aligned` ("on" by default) turns field-aligned
  interpolation and derivatives off for comparison; then parallel
  gradients are b·∇. For a plasma that doesn't vary in z the two give the
  same background to rounding (`test_field_aligned_cpu`).
  `hermes_snap_start` can count from the end, like `hermes_snap_end`.
- `python/check_hermes_bkg.py` checks the background as Flan built it,
  from an output file: the parallel ion flow into each target against
  the sound speed, and the direction of the ion grad-B drift relative to
  the X-point for each field direction, and prints where the impurities
  are by zone. Both checks were run on forward and reversed field runs,
  with and without field alignment.
- `python/plot_multizone.py` (formerly `plot_s3x_nz.py`, which still
  works) plots the impurity density or any background quantity of a
  SOLEDGE3X or HERMES run.
- [HERMES Backgrounds](../running/hermes.md) in Running Flan, and the
  options in [Input Options](../input-options.md#hermes).

**Later**: the GPU path, a time-averaged background, a sputtered source
at the targets, and field-aligned interpolation for turbulent SOLEDGE3X
runs. The BOUT++ Hermes-3 reader has its own plan, below.

## BOUT++ Hermes-3 reader

This is the plan for `hermes_format = "bout"`: reading the plasma from
the dump files a BOUT++ Hermes-3 run writes, instead of JAX-Hermes3
snapshots. The grid, the zones, the field-aligned interpolation, the
toroidal wedge and the conduction all stay as they are. Only how the
plasma gets from the files into the zoned arrays changes, plus two
physics points that JAX-Hermes3 runs don't raise (the flow in
electromagnetic runs, and the inductive parallel electric field).

### The dump files

- One file per MPI process, `BOUT.dmp.N.nc` for N = 0 to
  NXPE·NYPE - 1 (netCDF-4/HDF5). Each holds that process's block of
  MXSUB × MYSUB cells in (x, y) and all MZ z cells, **with** MXG x guard
  cells and MYG y guard cells on each side, for every output time:
  variables are (t, x, y, z).
- `PE_XIND` and `PE_YIND` give the block's place: global x of the block's
  first interior cell is MXG + PE_XIND·MXSUB, global y is PE_YIND·MYSUB.
  The global grid is NXPE·MXSUB + 2·MXG by NYPE·MYSUB, which must equal
  the grid file's `nx` (guards included) and `ny` (no y guards).
- To assemble one global field: from each block take the interior cells,
  plus its x guard cells on a side where the block is at the domain's x
  boundary (they are the grid's x guards, which the zones use). The y
  guards are always dropped, since the grid has none; at the targets they
  hold sheath boundary values, not plasma cells.
- `t_array` (in units of 1/Omega_ci = m_p / (e Bnorm)) holds the times.
  Every file has the same times. A restarted run starts its `t_array`
  where the previous run ended, and (with `append = false`) overwrites
  the earlier dumps.
- The normalization is in every file: `Nnorm`, `Tnorm`, `Bnorm`, `Cs0`,
  `Omega_ci`, `rho_s0`, and also `zperiod`, `MZ`, and the grid's
  separatrix indices. Each variable has `units` and `conversion`
  attributes. The ion mass and the conduction coefficients are not in the
  dumps; they are in `BOUT.settings` in the run directory.
- With `diagnose = true` the dumps also hold the derived quantities
  `Te`, `Ti`, `Vi`, `kappa_par_e`, `kappa_par_i`, and the conducted and
  total energy flows through the lower y faces (`efe_cond_ylow`,
  `efi_cond_ylow`, ...).
- Variables are stored one time per chunk. Opening a file costs much more
  than reading all of its times of one variable, so the reader has to
  open each file once and read every wanted time from it, not loop over
  times on the outside.

### The parallel flow in electromagnetic runs

With the `electromagnetic` component, Hermes-3 evolves the canonical
momentum, and the `NV` it writes is that: NVi = AA·n·V_i + Z·n·Apar for
each species (in normalized units, for electrons with Z = -1). So
NVi / (AA·n) isn't the flow there. The reader takes the flow from the
species' `V` (e.g., `Vi`) when it's in the dumps, and otherwise from
(NV - Z·n·Apar) / (AA·n), which is the same thing. It stops if neither
`V` nor `Apar` is there and the run is electromagnetic (the
`electromagnetic` component in `BOUT.settings`, or an `Apar` variable).

The parallel electric field in an electromagnetic run is
E_par = -grad_par phi - dApar/dt. Flan uses only the first term. With
`magnetic_flutter = false` the field lines (and so b and grad_par) are
the unperturbed ones from the grid, so nothing else changes.

### Phases

**Phase B0: conventions and data checks.** Done.
`python/check_hermes_conventions.py --format bout --run-dir <run>` runs
the checks of Phase 0 on a BOUT++ run, plus checks of the dumps
(section 6 of its output):

- It assembles global fields from the dump blocks as above, opening
  each file once (`--cache` keeps them in a `.npz`), and stops if a cell
  isn't filled exactly once or a file's times differ. Field-aligned
  differences across block boundaries are the same as inside blocks.
- It prints the `hermes_*` options to set from `BOUT.settings` (ion
  mass, conduction coefficients, flux limiter) and warns about settings
  Flan doesn't follow (ion charge other than 1, conduction collisions
  other than `multispecies`).
- The sources for the balance come from the dumped external sources
  (`Se_src`, and 3/2 `P_src` for each species) integrated over the
  volume, so no input expressions have to be parsed. Without recycling
  and before the density has settled, the particles reaching the
  targets can be several times the source; the balance then closes with
  the volume integral of `ddt(Ne)`.
- It fits NV - AA·N·V against N·Apar, confirming the canonical momentum
  in electromagnetic runs (slope Z to round-off).
- kappa_par from the script's Braginskii formulas (the same as Flan's
  `braginskii.cpp`) equals the dumped `kappa_par_e/i` to round-off. The
  conducted power through a y face, with the neighboring plane shifted
  spectrally onto the field line, matches `efe/efi_cond_ylow` upstream to
  well under 1%, and is correlated with it cell by cell everywhere. Next
  to the targets it is somewhat larger than the dumped power, for a
  reason still to be found (probably how Hermes-3 takes kappa or the
  face area next to the sheath). A linear shift in z, which the script
  used before, smooths the plane and misses the upstream power by
  several percent; the script now shifts spectrally, as BOUT++ and Flan
  do.
- It compares -dApar/dt (centered over neighboring output times, and
  over twice that to see whether the output interval resolves it) with
  -grad_par phi. In an electromagnetic turbulence run the inductive part
  can be a noticeable fraction of the electrostatic one in the
  fluctuations, and smaller in the toroidal average, so Phase B5 is
  worth doing as an option.

**Phase B1: split the snapshot reader, no change in behavior.** Done.
`read_snapshots` did three things: reads one snapshot, converts it to SI
in the zoned arrays, and finishes the background (time check, ghost
cells, E and the flows, heat fluxes). Split it so both formats share the
conversion and the finish:

- `Hermes::Frames`: the zoned ne, te, ti, upar and vp arrays and the
  times, with `store(t, x, y, k, Ne, Pe, Pi, V, phi)` taking normalized
  values at a global (x, y) cell (a lookup table from global (x, y) to
  (zone, i, j) built from `zone_ranges`). V is the parallel flow in units
  of c_s0 along +y, so each reader works out V its own way. `store`
  applies the normalization, `parallel_sign` and the Ne <= 0 rule, so
  they are in one place.
- `finish_background(frames, opts, bkg)`: everything after the loop in
  today's `read_snapshots`.
- The JAX reader becomes a loop over `.npz` files calling `store` with
  V = NVi / (AA Ne).
- The normalization (Nnorm, Tnorm, Bnorm, AA, zperiod) moves into a small
  `Hermes::Norm` struct filled per format (`norm_from_options` for
  JAX-Hermes3) and checked by `check_norm`, so the shared code doesn't
  read the options for it.

Check: `test_read_hermes_cpu` passes unchanged, and a JAX background read
before and after is the same bit for bit (a scratch dump driver, as for
the SOLEDGE3X refactor).

**Phase B2: the dump reader.** Done. `read_hermes_bout.cpp`:

- `read_bout_layout(dir)`: from `BOUT.dmp.0.nc`, NXPE, NYPE, MXSUB,
  MYSUB, MXG, MYG, MZ, `t_array`, Nnorm, Tnorm, Bnorm, `zperiod`, the
  separatrix indices and the variable names. Stops if `dir` has no
  `BOUT.dmp.0.nc`, NZPE > 1 (z split over processes isn't handled), the
  dumps' `nx`, `ny` don't follow from the layout, or a file is missing.
  `check_bout_layout` stops if the size or separatrix indices differ from
  the grid file's.
- `norm_from_bout`: Nnorm, Tnorm, Bnorm and `zperiod` from the dumps, the
  ion mass from `hermes_ion_mass_amu`. `hermes_nnorm`, `_tnorm`, `_bnorm`
  and `_zperiod` aren't used; a note is printed if one was changed from
  its default to a different value.
- Time selection: `hermes_snap_start`, `_end` and `_stride` index the
  dump times with the same rules as for `.npz` snapshots
  (`select_indices`, shared with `select_snapshots`). `hermes_snap_dir`
  is the run directory.
- `read_bout_dumps`: for each file, opened once, checks its place
  (`PE_XIND`, `PE_YIND`) and `t_array`, and reads every selected time of
  each variable as one strided hyperslab of its interior cells (and x
  guards at the domain edges), then `store`s them. It reads Ne, Pe,
  P<ion>, phi and the ion flow: V<ion> if dumped, else
  (NV<ion> - N Apar) / (AA N) in an electromagnetic run, else
  NV<ion> / (AA N). It stops if the files don't fill every global cell
  exactly once, and prints its progress.
- Before reading, the run directory, layout, normalization, the times
  picked and the memory of the plasma arrays are printed.
- Read on rank 0 and broadcast, as for every background. Reading files
  on several ranks and gathering is a later speed-up if needed.
- New option `hermes_ion_species` (default `"i"`), the ion species'
  name, giving the variable names `P<name>`, `NV<name>`, `V<name>`.
  Electrons are always `e`. `hermes_ion_mass_amu` and the
  `hermes_kappa_*` options are still needed, from `BOUT.settings`.

Check: the plasma Flan reads from a BOUT++ run (with and without a time
stride) equals the Python assembly in `check_hermes_conventions.py`
exactly for n, T and phi, and to the difference in proton mass between
CODATA versions (about 1e-8) for the flow. A JAX background is the same
bit for bit as before.

**Phase B3: unit tests.** Done. `HermesTest::write_bout_dumps`
(`tests/hermes_test_files.h`) writes the test plasma of the `.npz`
snapshots as dump files over 2 × 4 processes, with NaN in the y guard
cells and in the x guard cells between blocks, NV written as the
canonical momentum when Apar is written, and switches for each kind of
bad dump. Tests in `test_read_hermes_cpu`:

- `BoutMatchesJax`: the dumps and the snapshots give the same background
  bit for bit, with the normalization and `zperiod` from the dumps (the
  options are left at defaults that differ).
- `BoutFlowFromNV`: without V, the flow from (NV - N Apar) / (AA N) and,
  without Apar, from NV / (AA N), to rounding.
- `BoutTimeSelection`: strides and negative indices pick the right times.
- `BoutIgnoresNormOptions`, `BoutIonSpecies`.
- `BoutErrors`: a missing file, NZPE > 1, a file with other times, dumps
  of another grid, two files at the same place, a wrong ion species name,
  no flow variable, and times out of range each stop with a message.

Making the reader ignore Apar, or read the x guards between blocks, makes
these tests fail.

**Phase B4: validation on a BOUT++ run.**

- The plasma Flan reads, at a few times, against the Python assembler of
  Phase B0 (should be exact).
- Flan's kappa_par_e and kappa_par_i against the dumped `kappa_par_e`
  and `kappa_par_i`, and Flan's conducted heat flux times the y-face area
  against `efe_cond_ylow` and `efi_cond_ylow`. These are direct checks of
  `braginskii.cpp` and the field-aligned gradient on a run that wrote
  them, which the JAX snapshots couldn't give.
- `python/check_hermes_bkg.py` (target Mach numbers, grad-B drift
  direction) and a short tungsten run, as in Phase 6.

**Phase B5: inductive E_par (option, off by default).** Add
-dApar/dt to E_par. The derivative is taken when the dumps are read: for
each loaded time t_k, Apar is also read at the output times just before
and after it, and only the centered difference is kept. That works for
any `hermes_snap_stride`, needs no extra stored fields, and costs little
extra reading since the files are open anyway. A new option (e.g.
`hermes_inductive_epar = "off"/"on"`) turns it on; it stops if the run has
no Apar. Test: a synthetic Apar linear in time gives the expected E_par.

**Phase B6: docs.** `hermes_format`, `hermes_snap_dir`,
`hermes_snap_start` and the normalization options in
[Input Options](../input-options.md#hermes); a BOUT++ section in
[HERMES Backgrounds](../running/hermes.md) (which files, which options
come from `BOUT.settings`, memory and read time); this page's summary for
the HERMES-3 developers.

### Decisions to make

1. **Normalization from the dumps, not the options.** Recommended: the
   files are the authority, and taking them removes a way to get it wrong.
2. **Ion mass and conduction coefficients**: options (as now) or parsed
   from `BOUT.settings`. Recommended: options for now, with
   `check_hermes_conventions.py` printing the values from
   `BOUT.settings`, since the settings file is free-form text.
3. **Heat flux**: keep calculating q_par from Flan's Braginskii kappa
   (as for JAX), or add `hermes_heat_flux = "dumped_kappa"`, which uses
   the dumped `kappa_par_e/i` with Flan's grad_par T. Recommended: keep
   the calculation, use the dumped kappa to validate it (Phase B4), and
   add the option only if they disagree for a reason Flan can't fix.
4. **Inductive E_par**: Phase B0 found it large enough in an
   electromagnetic run to be worth having as an option (Phase B5).
5. **Several restart segments** (dumps in different directories, one
   after another in time): not in the first version. One directory per
   background; a list of directories could come later.

---

## Summary for the HERMES-3 developers

Flan is a Monte Carlo impurity transport code. It follows impurity ions
(mostly tungsten) through a fixed, time-dependent plasma background,
with the Lorentz force, Coulomb collisions with the main ions (including
the thermal force), and ionization and recombination from OpenADAS. The
impurities don't act back on the plasma. This summary describes how Flan
reads a HERMES background, so that the assumptions can be checked by the
people who know the code.

### What Flan reads

From the hypnotoad grid: `Rxy`, `Zxy` and the four cell corners, `Brxy`,
`Bzxy`, `Btxy`, `zShift`, `ShiftAngle`, `ixseps1`, `jyseps1_1`,
`jyseps2_2` and the metric (`J`, `g_22`, `dx`, `dy`). From each snapshot:
`Ne`, `Pe`, `Pi`, `NVi`, `phi` and `t`. From the run's input: Nnorm,
Tnorm, Bnorm, the ion `AA`, `zperiod` and the conduction coefficients.

### What Flan computes from them

- **Geometry**: each cell's (R, Z) corners and its toroidal extent give the
  cell in Cartesian coordinates, with (R, z, Z) taken as right handed.
- **Plasma**: n_e = Ne Nnorm, T_e = (Pe / Ne) Tnorm, T_i = (Pi / Ne) Tnorm
  (so n_i = n_e), V_par = NVi / (AA Ne) times c_s0 = sqrt(e Tnorm / m_p),
  and the potential phi Tnorm.
- **Fields**: B from `Brxy`, `Bzxy` and `Btxy`; E = -grad phi, with its
  parallel part -grad_par phi taken along the field line as below.
- **Ion flow**: V_par b plus the E x B and grad-B drifts.
- **Parallel heat fluxes** (for the thermal force), as in
  `braginskii_conduction`: q_par = -kappa grad_par T with
  kappa = c P tau / m and tau = 1 / Σν, where Σν is the sum of the
  species' Braginskii collision frequencies (e-e and e-i for electrons,
  i-i and i-e for ions, as in `braginskii_collisions`). The flux limiter
  (`kappa_limit_alpha`) is off unless it's set.
- **Parallel derivatives and interpolation in y** follow the field line:
  going from plane y to y + 1, z changes by Δ`zShift`, plus `ShiftAngle`
  across the core branch cut. For derivatives the neighboring planes are
  shifted in z with a Fourier transform, and
  grad_par f = (B·∇y / B) ∂f/∂y in field-aligned coordinates, as in
  BOUT++'s `Grad_par`. Interpolation at the particles is linear in z.
- **Time**: Flan interpolates linearly in time between the loaded
  snapshots.

### Assumptions

1. z is the lab-frame toroidal angle (shifted metric), and the field-line
   shift between planes is Δ`zShift` as above. On the data this lines up
   the turbulence between neighboring planes everywhere except right next
   to the X-point.
2. (R, z, Z) is right handed. The grid's `J` and `bxcvx` agree with this,
   and with it the ion grad-B drift in the "forward" grid points toward the
   X-point.
3. NV is the momentum density AA·n·V, with V positive toward +y, and +y
   along the poloidal part of b (the poloidal field points along +y).
4. One ion species of charge 1, quasineutral, so n_i = n_e.
5. The domain covers 1/`zperiod` of the torus and is periodic in z. A
   particle crossing the z boundary reappears rotated by 2π/`zperiod`, and
   its velocity is rotated with it.
6. The y ends of the legs are material targets. The outer SOL and PFR x
   boundaries absorb impurities, and the core boundary is the edge of the
   simulated region.
7. Without neutrals in the run, Flan sees only the ions and electrons.
   Impurity ionization and recombination use n_e and T_e.
8. Snapshot fields are cell-centered values with the x guard cells
   included and no y guard cells.

### Questions for the HERMES-3 developers

1. Does JAX-Hermes3's parallel conduction match Hermes-3's
   `braginskii_conduction` with the default `multispecies` collisions and
   no flux limiter? If not, what should Flan use, or could q_par be
   written to the snapshots?
2. Are assumptions 1 to 3 right for both field directions, including the
   sign of `ShiftAngle` at the core branch cut?
3. Some of the particle and energy sources don't reach the targets
   through the parallel fluxes alone. Is the rest the radial flux through
   the outer SOL and PFR boundaries, or something else that Flan should
   include (e.g., a sheath heat flux different from the bulk fluxes)?
4. For BOUT++ Hermes-3 output: is the plan in
   [BOUT++ Hermes-3 reader](#bout-hermes-3-reader) right about the guard
   cells, the canonical momentum written as `NV` in electromagnetic runs,
   and `V` being the flow? Is dApar/dt usually small enough to leave out
   of E_par for impurities?
5. Does JAX-Hermes3 shift between planes with FFTs (as BOUT++'s
   `ShiftedMetric` does) and take grad_par as above? If it interpolates
   differently, Flan's E∥ and q_par would differ slightly from the run's.
