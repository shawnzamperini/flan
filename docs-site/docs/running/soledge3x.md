# SOLEDGE3X Backgrounds

Flan can follow impurities in a SOLEDGE3X plasma background, set with
`bkg_source = "soledge3x"`. This page explains how SOLEDGE3X backgrounds
work in Flan and what is different from Gkeyll backgrounds. An example
input is `regression/s3x_example.cpp`, and the options are listed under
[SOLEDGE3X](../input-options.md#soledge3x) on the Input Options page.

!!! note "GPUs"
    SOLEDGE3X backgrounds run on GPUs (`use_gpu = "cuda"`) as well as the
    CPU, except for the thermal force (with heat fluxes in the plasma
    files, set `imp_thermal_force = "off"`), the sputtered source and
    self-sputtering, which stop, and the sheath field, which a GPU run
    goes without. These are CPU only so far.

## Files

Everything is read from the directory `s3x_dir`:

- The mesh (`s3x_mesh_file`), with the cell corners, magnetic field and
  wall mask of every zone, and how the zones connect. Only axisymmetric
  meshes are supported.
- The reference file (`s3x_ref_file`), with the values SOLEDGE3X
  normalizes its output by.
- The plasma files `plasma_%05d.h5` from `s3x_frame_start` to
  `s3x_frame_end`, with the electron and ion density and temperature, the
  potential and the parallel ion flux of each frame.

The plasma is unnormalized with the reference values: densities by n0,
temperatures by T0, the potential by phi0, the parallel ion flow is G / n
times c0, and the time is tsimu times tau0. SOLEDGE3X's reference mass is
1 amu, so c0 and tau0 are not the main ion's sound speed and inverse
cyclotron frequency. The electric field is minus the gradient of the
potential, and the background ion flow is the parallel flow plus the
E x B and grad-B drifts, the same as for Gkeyll backgrounds.

If the plasma files have the fluxes on the cell faces (fluxn and fluxE,
for both species and in every frame), Flan also computes the parallel
conductive heat fluxes of the ions and electrons from the fluxes on the
phi faces, and uses them for the thermal force in the collisions (see
[Thermal Force](../physics/thermal-force.md) and `imp_thermal_force`). The
constant that converts the saved fluxes to physical units is fitted
against the parallel flow and printed for each species at the start of
the run; it must be within 5% of -2 pi c0 / a0. The heat fluxes are saved
with the background as `q_par_i` and `q_par_e` (W/m²). Without the face
fluxes, a warning says the thermal force isn't included.

SOLEDGE3X's (R, Z, phi) is right-handed, so phi increases clockwise seen
from above.

## Zones

SOLEDGE3X splits its grid into zones to handle the X-point, e.g., the
core, the scrape-off layer, the private flux region and the divertor
legs. Flan keeps the same zones, numbered from 1 like SOLEDGE3X, and every
particle knows which zone it is in. Within a zone, the computational
coordinates are:

- x: the psi cell index, from 0 to the zone's number of psi cells.
- y: the theta cell index, from 0 to the zone's number of theta cells.
- z: the toroidal angle phi, in radians, between the zone's
  `phicorners`. If they span less than the whole torus, the run is a
  toroidal wedge: a particle that leaves through one end comes back in
  at the other, with its velocity rotated by the wedge's width.

So x = 5.5, for instance, is the middle of the sixth psi cell of the zone.
The starting position of particles (`imp_xstart_opt`, etc.) is given in
these coordinates, in the zone set by `imp_zone_start_val`.

## Boundaries

What happens at each face of each zone comes from SOLEDGE3X's mesh, not
from the bound type options (`min_xbound_type`, etc.), which are only for
Gkeyll and slab backgrounds. A particle that crosses a face shared with
another zone continues in that zone. Faces at a divertor target are wall
faces, and faces at the outer psi boundary are absorbing faces. Both
absorb the particle, but they are counted separately in the particle
fates printed at the end of a run, and wall faces are material surfaces
(e.g., for [sputtering](../physics/sputtering.md)).

SOLEDGE3X's neighbor codes map to Flan's face codes as:

| SOLEDGE3X | Meaning | Flan |
|---|---|---|
| -1 | Target | -4 (wall) |
| -2 | Core | -2 (core) |
| -3 | Outer psi boundary | -1 (absorbing) |

Flan uses -4 rather than -3 for wall faces so that -3 never means
different things in the mesh file and in Flan's output.

A particle that crosses into the core is put back at a random place on
the core boundary, picked between the core zones in proportion to their
number of cells. With `s3x_chi_absorb = "on"`, a particle that ends a time
step in a cell inside the wall is absorbed.

The buffers (`imp_xbound_buffer`, etc.) still apply at absorbing, wall and
core faces, but in cell-index units for x and y (e.g., 0.5 is half a
cell).

## Ghost cells

Fields are interpolated across the faces shared by zones using ghost
cells, two cells deep on every side, filled from the neighboring zones.
SOLEDGE3X's plasma files have ghost cells of their own, and when the
background is loaded, Flan prints how they compare with Flan's. They
should match exactly everywhere except at the boundaries, where SOLEDGE3X
stores its boundary conditions, and at the X-point corners, which
SOLEDGE3X computes in its own way. Flan uses the average of the cells on
either side there. Only interpolation in the cells touching the X-point
is affected.

## Output

With more than one zone, the netCDF file is laid out differently from a
Gkeyll one, since each zone has its own grid:

- `input` and `geometry/time` are the same as usual.
- `zones` has the number of cells in each zone (`nx`, `ny`, `nz`) and
  `faces`, what is across each face of each zone: a zone number (from 1),
  -1 for an absorbing face, -2 for a core face or -4 for a wall face. The
  rows are the xmin, xmax, ymin, ymax, zmin and zmax faces.
- `zone1`, `zone2`, ... each have their own x, y, z and grid dimensions,
  and `geometry`, `background` and `output` groups with the usual
  variables for that zone's cells. The geometry also has the wall mask
  `chi` and `R_corners` and `Z_corners`, the (R, Z) cell corners in
  meters, for plotting.
- `sputtering`, with `imp_source_opt = "sputter"`, has the wall surfaces
  and the rate the background sputters atoms off them. See
  [Sputtering](../physics/sputtering.md).

The output is large for a full SOLEDGE3X mesh, since every zone stores
its background and output for every frame. It grows with the number of
cells and frames.
