=========================================================================================
SOLEDGE3X Backgrounds
=========================================================================================

Flan can follow impurities in a SOLEDGE3X plasma background, set with :ref:`bkg_source <bkg_source>` = "soledge3x". This page explains how SOLEDGE3X backgrounds work in Flan and what is different from Gkeyll backgrounds. An example input is ``regression/s3x_example.cpp``.

Only the CPU version of Flan can run SOLEDGE3X backgrounds so far (``use_gpu`` = "off").

Files
-----

Everything is read from :ref:`s3x_dir <s3x_dir>`:

  - The mesh (:ref:`s3x_mesh_file <s3x_mesh_file>`), with the cell corners, magnetic field and wall mask of every zone, and how the zones connect. Only axisymmetric meshes are supported.
  - The reference file (:ref:`s3x_ref_file <s3x_ref_file>`), with the values SOLEDGE3X normalizes its output by.
  - The plasma files ``plasma_%05d.h5`` from :ref:`s3x_frame_start <s3x_frame_start>` to :ref:`s3x_frame_end <s3x_frame_end>`, with the electron and ion density and temperature, the potential and the parallel ion flux of each frame.

The plasma is unnormalized with the reference values: densities by n0, temperatures by T0, the potential by phi0, the parallel ion flow is G / n times c0, and the time is tsimu times tau0. SOLEDGE3X's reference mass is 1 amu, so c0 and tau0 are not the main ion's sound speed and inverse cyclotron frequency. The electric field is minus the gradient of the potential, and the background ion flow is the parallel flow plus the E x B and grad B drifts, the same as for Gkeyll backgrounds.

SOLEDGE3X's (R, Z, phi) is right-handed, so phi increases clockwise seen from above.

Zones
-----

SOLEDGE3X splits its grid into zones to handle the X-point, e.g., the core, the scrape-off layer, the private flux region and the divertor legs. Flan keeps the same zones, numbered from 1 like SOLEDGE3X, and every particle knows which zone it is in. Within a zone, the computational coordinates are:

  - x: the psi cell index, from 0 to the zone's number of psi cells.
  - y: the theta cell index, from 0 to the zone's number of theta cells.
  - z: the toroidal angle phi, in radians.

So x = 5.5, for instance, is the middle of the sixth psi cell of the zone. The starting position of particles (:ref:`imp_xstart_opt <imp_xstart_opt>`, etc.) is given in these coordinates, in the zone set by :ref:`imp_zone_start_val <imp_zone_start_val>`.

Boundaries
----------

What happens at each face of each zone comes from SOLEDGE3X's mesh, not from the bound_type options (:ref:`min_xbound_type <min_xbound_type>`, etc.), which are only for Gkeyll backgrounds. A particle that crosses a face shared with another zone continues in that zone. Faces at a divertor target are wall faces, and faces at the outer psi boundary are absorbing faces. Both absorb the particle, but they are counted separately in the particle fates printed at the end of a run, and wall faces are material surfaces (e.g., for sputtering). SOLEDGE3X's neighbor codes map to Flan's face codes as -1 (target) to -4 (wall), -2 (core) to -2 (core) and -3 (outer psi boundary) to -1 (absorbing). Flan uses -4 rather than -3 for wall faces so that -3 never means different things in the mesh file and in Flan's output. A particle that crosses into the core is put back at a random place on the core boundary, picked between the core zones in proportion to their number of cells. With :ref:`s3x_chi_absorb <s3x_chi_absorb>` = "on", a particle that ends a time step in a cell inside the wall is absorbed.

The buffers (:ref:`imp_xbound_buffer <imp_xbound_buffer>`, etc.) still apply at absorbing, wall and core faces, but in cell-index units for x and y (e.g., 0.5 is half a cell).

Ghost Cells
-----------

Fields are interpolated across the faces shared by zones using ghost cells, two cells deep on every side, filled from the neighboring zones. SOLEDGE3X's plasma files have ghost cells of their own, and when the background is loaded, Flan prints how they compare with Flan's. They should match exactly everywhere except at the boundaries, where SOLEDGE3X stores its boundary conditions, and at the X-point corners, which SOLEDGE3X computes in its own way. Flan uses the average of the cells on either side there. Only interpolation in the cells touching the X-point is affected.

Output
------

With more than one zone, the netCDF file is laid out differently from a Gkeyll one, since each zone has its own grid:

  - ``input`` and ``geometry/time`` are the same as usual.
  - ``zones`` has the number of cells in each zone (nx, ny, nz) and faces, what is across each face of each zone: a zone number (from 1), -1 for an absorbing face, -2 for a core face or -4 for a wall face. The rows are the xmin, xmax, ymin, ymax, zmin and zmax faces.
  - ``zone1``, ``zone2``, ... each have their own x, y, z and grid dimensions, and geometry, background and output groups with the usual variables for that zone's cells. The geometry also has the wall mask chi and R_corners and Z_corners, the (R, Z) cell corners in meters, for plotting.
  - ``sputtering``, with :ref:`imp_source_opt <imp_source_opt>` = "sputter", has the wall surfaces and the rate the background sputters atoms off them. See the Sputtering page of the documentation website.

The output is large for a full SOLEDGE3X mesh, about 1.7 GB for two frames of the DIII-D example.
