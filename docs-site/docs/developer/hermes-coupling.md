# HERMES Coupling Plan

This page is the development plan for reading HERMES plasma backgrounds
into Flan, the same way SOLEDGE3X backgrounds are read
(`bkg_source = "hermes"`). The first version reads output from
JAX-Hermes3, a Python/JAX reimplementation of Hermes-3. A reader for
BOUT++ Hermes-3 dump files is planned but not started. The last section
is a short summary for the HERMES-3 developers of how the coupling works
and what it assumes.

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
(`BOUT.dmp.N.nc`) holding every time, with the normalization in the file
and the variables named per species (`Nd+`, `Pd+`, `NVd+`, ...). The grid
file and the physics are the same, so only the snapshot reader differs.

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

**Phase 3: field-aligned interpolation and gradients.** On a shifted-metric
grid a field line moves many toroidal cells between neighboring y planes,
so interpolating or differentiating in y at fixed z mixes unrelated
turbulent structures. Each zone gets an optional `zShift(x, y)` table.
With it, `Utilities::interp_vec4d` takes the y + 1 plane at
z + Δ`zShift` (interpolated in z, periodic), and y derivatives (for
E = -grad phi and grad_par T) are taken the same way. Without the table
nothing changes. Unit test: a field-aligned analytic function must
interpolate and differentiate exactly.

**Phase 4: toroidal period.** Each zone gets an optional toroidal period
(2π / `zperiod` for HERMES, the span of `phicorners` for SOLEDGE3X). When
a particle wraps in z on a zone with a period under 2π, its Cartesian
velocity (and position, if stored) is rotated by the period, since the
point it wraps to is rotated by that angle. Until then, the SOLEDGE3X
reader stops if `phicorners` doesn't span 2π. Tested in
`test_boundary_cpu` and with a SOLEDGE3X wedge test file.

**Phase 5: thermal force.** q_par = -kappa_par grad_par T for each
species, with kappa_par as in Hermes-3's `braginskii_conduction` and the
field-aligned derivative of Phase 3. The heat fluxes are saved with the
background like the SOLEDGE3X ones.

**Phase 6: validation and output.** `regression/hermes_example.cpp` like
`s3x_example.cpp`, starting from a point source; a check that V/cs is
about 1 or more at the targets and that the grad-B drift has the right
direction for each field direction; agreement with and without field
alignment for fields that don't vary in z; plotting through a
generalized `plot_s3x_nz.py`; and a page in Running Flan plus the options
in [Input Options](../input-options.md).

**Later**: the GPU path, a time-averaged background, the BOUT++ Hermes-3
reader, a sputtered source at the targets, and field-aligned
interpolation for turbulent SOLEDGE3X runs.

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
- **Fields**: B from `Brxy`, `Bzxy` and `Btxy`; E = -grad phi.
- **Ion flow**: V_par b plus the E x B and grad-B drifts.
- **Parallel heat fluxes** (for the thermal force), as in
  `braginskii_conduction`: q_par = -kappa grad_par T with
  kappa = c P tau / m and tau = 1 / Σν, where Σν is the sum of the
  species' Braginskii collision frequencies (e-e and e-i for electrons,
  i-i and i-e for ions, as in `braginskii_collisions`). No flux limiter.
- **Parallel derivatives and interpolation in y** follow the field line:
  going from plane y to y + 1, z changes by Δ`zShift`, plus `ShiftAngle`
  across the core branch cut.
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
4. For BOUT++ Hermes-3 output: are there any differences in variable
   names, normalization or guard cells from the above that a reader would
   need to handle?
