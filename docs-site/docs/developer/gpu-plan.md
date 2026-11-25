# GPU Plan: Multi-Zone Backgrounds, Sputtering and the Sheath

This page is the development plan for bringing the GPU (CUDA) path up to
the CPU path. The SOLEDGE3X and HERMES backgrounds, the toroidal wedge,
field-aligned interpolation, the thermal force, the sputtered source, the
sheath and self-sputtering were all written for the CPU first. The work is
on the branch `gpu-bkgs-sheath`.

## Why

The GPU path was written when Flan only had Gkeyll backgrounds: one
logically rectangular zone, particles started from the `imp_*start_opt`
options, and absorbing, periodic or core boundaries. Everything added
since then stops a GPU run, or is quietly left out of it. The runs that
need the most particles (sputtered sources on SOLEDGE3X walls and HERMES
targets, with the sheath and self-sputtering) are the ones that can't use
the GPUs.

## What the GPU can't do now

Each item below is either stopped with an error or left out of a GPU
run today:

| Feature | On a GPU run now | Where |
|---|---|---|
| Sputtered source (`imp_source_opt = "sputter"`) | Runs (since G5) | `sputtering_hd.h`, `cuda/slots.cu` |
| Deposition and prompt tallies, impact histograms | Run (since G5) | `sputtering_hd.h`, `cuda/boundary.cu` |
| Sheath field and near-wall sub-steps | Runs (since G6) | `sheath_hd.h`, `cuda/boris.cu` |
| Self-sputtering (`sput_self`) | Stops (since G5) | `flan.cpp`, `Sputtering::SelfSource::rebirth` |
| `imp_iz_max_charge` (test option) | Runs (since G5) | `cuda/openadas.cu` |
| Boundary fates: `self_sputtered` | 0 (no self-sputtering) | `follow_impurities` |
| `Slots` fields `prompt_by`, `generation` | On the device since G5 (`zone` since G1) | `slots_device.h` |

Things that are host-side setup and need nothing on the GPU: reading the
backgrounds (Gkeyll, SOLEDGE3X, both HERMES formats), the field-aligned
derivatives, the conductive heat fluxes and the inductive E_par (all
computed into the background before the run), `Sputtering::Surface`, the
source's rates and the D+ impact table (`Sheath::ImpactTable`, built at
startup), and saving the results.

## What the GPU has now

- Device structs of raw pointers for the background, slots, statistics,
  OpenADAS and options (see [Code Structure](code-structure.md)). Before
  G1 the background's 1-D grids (times, cell centers and edges) were in
  constant memory, `MAX_1D = 500` each, so only one background fit on a
  GPU; G1 moved them to global memory with the rest of the zones.
- Before G1 the background data on the device was a single zone with no
  ghost cells, which has the same layout as a `Vector4D`. A `ZonedVector`
  is already flat (zone after zone, with per-zone offsets and strides),
  so G1 copies the multi-zone data as it is.
- Interpolation through an 8- or 16-point stencil built once per particle
  per step (`include/cuda/interpolate.cuh`), clamped at the grid's edges.
- Kernels for every main loop step, with the Nanbu collisions (including
  sub-cycling), OpenADAS ionization and recombination, the Boris push in
  curvilinear coordinates, absorbing/periodic/core boundaries, filling
  slots from the `imp_*start_opt` options, and the statistics.
- Since G2, every kernel works in any zone, with the plasma interpolated
  along the field line in zones with a zShift table, and births in
  `imp_zone_start_val`. Since G3 the boundaries are the CPU's own
  per-particle routines (`boundary_hd.h`), so SOLEDGE3X and HERMES runs,
  toroidal wedges and the `lcfs` boundary run on the GPU, and the fates are
  counted there.
- Since G4 the thermal force (heat-flux-skewed sampling) runs on the GPU
  with the CPU's sampler (`collisions_hd.h`).
- Since G5 the sputtered source, the charge limit, the prompt window and
  the deposition tallies run on the GPU with the CPU's routines
  (`sputtering_hd.h`).
- Since G6 the sheath field and the near-wall sub-steps run on the GPU
  with the CPU's routines (`sheath_hd.h`).
- One host thread per GPU, slots split between ranks and then GPUs.
- Tests `*_cuda`, and the regression variants `z10_gpu` (friction and
  thermal force), `gradb_gpu` (drifts) and `p01_gpu`, `p1_gpu`, and
  any sheath case followed by `_gpu` (prompt redeposition).

## Approach

### One code path for one zone and many

The GPU background becomes a device copy of the zones and the
`ZonedVector`s: per-zone cell counts, ghost width, offsets and strides,
the centers and edges of every zone in global memory, the face codes,
the toroidal, Cartesian and wedge flags, and the zShift tables. A Gkeyll
or slab background is one zone of this, with no special case. Constant
memory goes away (except perhaps a few scalars): with particles spread
over many cells, the reads of a warp aren't uniform, and constant memory
serializes those. Read-only global memory (`const __restrict__`, `__ldg`)
should be as fast. G1 checks that the single-zone runs don't get slower.

### Share per-particle code between the CPU and the GPU

Most GPU routines today are copies of the CPU ones (`nanbu_s_a.h` and
`nanbu_s_a.cuh` are the same 3357 lines). CLAUDE.md asks that the two
stay in step, and copies drift apart. The new per-particle code (the
sheath profile and field, the sub-step count, the skew sampler, zone
crossings, the wedge rotation, the Thompson and cosine sampling, the
yield table lookup) has no state of its own, so it can be written once as
`FLAN_HD` functions (`__host__ __device__` with CUDA, nothing without) in
headers that take plain pointers and sizes, and be called by both the CPU
code and the kernels. The CPU classes keep their interfaces and call the
shared functions inside. The existing duplicated code is only moved over
where a phase touches it anyway.

### Everything stays on the device during the loop

The main loop has no host round trips except the `all_dead` check. That
stays true: the sputtered source samples births on the device, deposition
and impact tallies are device arrays added to with `atomicAdd` (double
precision `atomicAdd` is native from sm_60, and the A100s are sm_80),
copied back and added to the host objects after the threads join, like
the statistics are now.

### Checking

Every phase is checked three ways:

1. **Deterministic tests** (`*_cuda`): the shared function or kernel
   against the CPU version on the same inputs (interpolation at given
   points, find-cell, a list of zone crossings, the sheath field, sub-step
   counts, a birth from given random numbers). These should agree to
   rounding.
2. **Regression variants** (`*_gpu`): each slab regression case touched
   by a phase gets a GPU variant, checked by the same checker script.
   CPU and GPU results won't be identical (different random number
   streams), so they agree within the Monte Carlo error.
3. **End-to-end runs**: the SOLEDGE3X and HERMES example cases on both,
   comparing density, deposition and fates within the Monte Carlo error,
   and the run time.

## Phases

**G0. Baseline.** Done.

- `build_gpu` from the branch: every `*_cuda` test passes, and so do
  `z10_gpu` and `gradb_gpu` (`check_collisions.py`, `check_drifts.py`).
  `MPITest` needs `srun`, so it was not run.
- The `FLAN_HD` macro (`include/host_device.h`, on `__CUDACC__`, not
  `USE_CUDA`, since a GPU build's `.cpp` files go through the host
  compiler) and the first shared function, `Utilities::trilinear_box`
  (`include/utilities_hd.h`). The CPU's `trilinear_interpolate` calls it
  with `x1 - x0`, and the GPU stencils call it with their signed box sizes,
  so neither side's arithmetic changed: `gradb` and `gradb_gpu` give the
  same output bit for bit as before. The GPU copy in
  `cuda/utilities.cuh` is gone. A NaN weight (a box of size 0) is now
  clamped to 0 on both sides; before, the CPU kept the NaN.
  `test_trilinear_cuda` compares the device with the host on random boxes
  of both signs, points outside the box and boxes of size 0.
- The device and host don't agree bit for bit, only to rounding: nvcc
  fuses `a + t * (b - a)` into a fused multiply-add. So CPU against GPU is
  always a comparison to rounding (deterministic tests) or within the
  Monte Carlo error (runs), as decision 7 says.
- GPU runs weren't reproducible, even with a fixed `seed`.
  `compact_dead_slots` listed the dead slots in an order set by atomics,
  and when there are more dead slots than particles left to launch, that
  order decides which slots (each with its own random number generator)
  get them. Replaced by `cub::DeviceSelect::If` over the slot indices,
  which lists them in increasing order; its work space is kept in
  `SlotsDevice::select_temp`. Two `z10_gpu` runs with a fixed seed now
  have the same counts, and the velocity sums differ by about 1e-11 m/s
  (the order of the statistics' `atomicAdd`s). This is what makes
  G1 checkable: a refactor that doesn't change the arithmetic must give
  the same counts and the same sums to rounding. The regression cases
  don't set `seed`, so these checks use copies with one.
- Timing on the login node's A100 is not a baseline: other users' jobs
  kept it at 100%, and `gradb_gpu` (one particle, 10⁴ steps) took 88 s,
  nearly all in the two host syncs per step (fill slots and the all-dead
  check), against 3 s on 16 CPU threads. `z10_gpu` (collisions) took
  9.3 s against 7.6 minutes on 16 CPU threads. Single-zone timings for
  G1 must come from a GPU compute node.
- Baseline on a GPU node (4 A100-SXM4-80GB, one rank, 64 OpenMP
  threads for the CPU runs, every GPU run on all 4 GPUs), the main loop
  ("Following impurities") with a fixed seed:

    | Case | CPU | GPU |
    |---|---|---|
    | `gradb` (1 particle, 10⁴ steps) | 2.08 s | 2.61 s |
    | `z10` (friction, collisions) | 63 s | 2.83 s |
    | Gkeyll, W over the whole volume, collisions and ionization, 50 frames, 1e5 particles | 90.6 s | 3.74 s |
    | Same, 1e6 particles | — | 10.3 s |

  The Gkeyll CPU and GPU runs agree (summed density 1812.0 and 1811.7,
  1815.5 with 1e6). Reading the background took 11.7 minutes in the
  first run, which converts the Gkeyll files into
  `bkg_from_pgkyl_*.bin`, and 1 s in the later ones, which load that. G1 is checked against these totals with the same job
  (`$SCRATCH/flan_gpu/g0_node/timing.sbatch`).
- Two things about timing GPU runs, for G8: the per-section times in the
  timer summary are launch times, since the launches are asynchronous,
  and the kernels' time shows up at the first sync ("Filling slots",
  "Checking particle state"); only the total is a real time. And with
  several GPUs the four threads add into one timer, so the sections sum
  to more than the total. Setting `CUDA_VISIBLE_DEVICES` before `srun`
  didn't limit a run to one GPU (`srun --gpus=1` would).
- `mkdocs.yml`: `FLAN_HD=` in Doxygen's `PREDEFINED` (mkdocs isn't
  installed here, so not built).

**G1. Zones on the device.** Done.

- One struct for the host and the device instead of a separate
  `ZonesDevice`: `Zones::ZonesView` (`include/zones_hd.h`) holds per-zone
  nx, ny, nz, the ghost width, the face codes (`num_faces` per zone), the
  centers (with the ghost cells' centers, NaN past a boundary face) and
  edges of every zone in flat arrays with per-zone offsets, the
  `ZonedVector` zone offsets (one set per background, since every field
  has the same layout), the toroidal, Cartesian and wedge flags, and the
  zShift tables. `Zones` keeps a host copy (`FlatZones`, rebuilt by every
  function that changes the zones; copying a `Zones` points the copy's
  view at its own arrays), `Zones::view` returns it, and
  `Zones::to_device`/`free_zones` copy it to a GPU and free it. The face
  codes and `Face` moved to `zones_hd.h` with it.
- Shared `FLAN_HD` functions on a view: `cell_index`, `neighbor_index`,
  `bracket`, `center` and `flat_index`. The CPU's `Zones::neighbor_index`
  and `Zones::bracket` call them, and `Zones::center` reads the same
  table within the ghost width (and steps across faces past it, as
  before). The table is filled by the old stepping code, so the CPU's
  values didn't change.
- `BackgroundDevice` holds the view, the times and their bounds in global
  memory, and the fields as the flat `ZonedVector` data, ghosts included.
  `to_device` checks every field it copies has the size the zones and
  times give. `device_constants.cu`/`.cuh` and `MAX_1D` are gone. The
  multi-zone and wedge stops moved from `to_device` to `flan.cpp`, so
  tests can copy any background; they come out in G3.
- `SlotsDevice::zone`, copied both ways. GPU births go in zone 0 until
  G2 (`imp_zone_start_val`), with the range and cell from the zone's
  edges. Every kernel takes its zone from the slot: the stencils
  (`build_stencil_3d`/`4d` take a zone, and store a zone/frame base plus
  int offsets per corner instead of 16 flat ints), the cell finding, the
  births, and the (still single-zone) boundary kernel's bounds.
- Tests: `test_zones_cuda` (cells, neighbors, brackets and centers on the
  device against the host for 20000 points on three zones with uneven
  cells, a neighbor thinner than the ghost width, periodic, core and wall
  faces; ghost centers against `Zones::center`; `flat_index` against
  `ZonedVector::calc_index`), and `test_interpolate_cuda` rewritten: the
  slab case, and two zones with ghost cells where 4000 points interpolate
  through the shared face exactly and agree with the CPU's `interp_ne`.
  The SOLEDGE3X and HERMES test grids are written with HDF5, which the
  `.cu` tests can't include, so the zones are built in the test.
- Single zone unchanged: with a fixed seed, `z10_gpu` has the same counts
  as G0 and its velocity sums differ by about 1e-11 m/s (as two G0 runs
  do); `gradb_gpu` and `gradb` are bit for bit the same as G0. CPU and GPU
  `ctest` pass (except `MPITest`, which needs `srun`).
- Timing: the G0 job again, `$SCRATCH/flan_gpu/g1_node/timing.sbatch`
  (same cases and configuration, plus one Gkeyll 1e6 run with
  `srun --gpus=1`). This node had the 40 GB A100s, G0's the 80 GB ones.
  Main loop, G0 → G1:

    | Case | CPU | GPU (4 GPUs) |
    |---|---|---|
    | `gradb` | 2.08 → 2.47 s | 2.61 → 2.48 s |
    | `z10` | 63 → 58 s | 2.83 → 2.75 s |
    | Gkeyll, 1e5 particles | 90.6 → 77.4 s | 3.74 → 3.76 s |
    | Gkeyll, 1e6 particles | — | 10.3 → 11.0 s |

  The Gkeyll 1e6 run on one GPU took 17.5 s. Only the 1e6 run is slower,
  by about 6% (G0's two 1e6 runs differed by 1%); the CPU column moves by
  more than that between jobs. That isn't enough to bring back a
  single-zone fast path (decision 1); G8 looks at it again. The seeded
  outputs match G0's: `gradb` and `gradb_gpu` bit for bit, `z10`, `z10_gpu`
  and the Gkeyll 1e5 CPU and GPU runs with the same counts and velocity
  sums within 1e-11 m/s. The Gkeyll 1e6 runs on one and four GPUs write
  the same file, so the 4-GPU 1e6 output wasn't compared.

**G2. Multi-zone push, collisions, ionization and statistics.** Done.

- The stencils were already per zone with ghost cells at shared faces
  (G1). What G2 adds is field-aligned interpolation on the device:
  `Utilities::interp_aligned_flat` (`include/interp_hd.h`, `FLAN_HD`) does
  the work of `Utilities::interp_aligned` on the fields' flat data,
  scalars and vectors (each sample rotated from its cell's toroidal angle
  to the point's), and the CPU's `interp_aligned` now calls it. It needs
  `has_z_shift`, `z_shift` and `wrap_z` on a `ZonesView`, which moved to
  `zones_hd.h` (the CPU's `Zones::z_shift` and `wrap_z` call them), and
  `bracket_sorted` for the times. `BackgroundDevice` has `phi_sign`.
- `Interpolate::interp_plasma` and `interp_plasma_vec3`
  (`cuda/interpolate.cuh`): the shared aligned interpolation in a zone with
  a zShift table, the stencil elsewhere. Used for ne, Te, Ti and the flow in
  the collision kernel, ne and Te in ionization/recombination, Ti in the
  births (`imp_temp_start_opt = "main_ion"`), and E and |E| in the Boris
  kernel, which rescales the aligned E to the aligned |E| as
  `push_one_cpu` does. B and the basis vectors keep the stencil, as on the
  CPU. (`interp_tang` isn't used by either push, so there's nothing to
  port.)
- `StatisticsDevice` is per zone: the flat `ZonedVector` data (no ghost
  cells) with each zone's offset and nx, ny, nz on the device, indexed with
  `ImpurityStats::index`. `to_device` and `add_stats_device` no longer go
  through the single-zone `get_dim2`.
- Births start in `imp_zone_start_val` (`OptionsDevice::zone_start`), and
  `fill_slots_kernel` writes the birth's time frame and cell into the slot
  as the CPU does, instead of zeros (the main loop finds them again right
  after anyway).
- The multi-zone stop in `flan.cpp` stays until G3: every kernel but the
  boundary check now handles zones, and a particle leaving its zone needs
  the boundary routine.
- Tests, all on the HERMES test grid (six zones with zShift, toroidal) and
  the SOLEDGE3X test mesh (three zones, no zShift), with particles spread
  over every zone (`tests/multizone_cases.h`):
  `test_boris_cuda` (`update_velocity_gpu` against `update_velocity_cpu`,
  Cartesian and computational velocities to 1e-10 of their largest
  values; with the aligned E turned off on the GPU the HERMES case is off
  by 1e-5 to 1e-4), `test_plasma_interp_cuda` (ne, Te, Ti, |E| and the flow
  through `interp_plasma` against `Background::interp_ne` etc., to
  1e-11; its kernel is in `plasma_interp_kernel.cu`, so the test itself is
  compiled by the host compiler and can write the test files), and
  `test_stats_cuda` (`record_stats_gpu` against `record_stats_cpu`, equal
  counts and sums to 1e-12; births in zone 3 of the SOLEDGE3X mesh stay
  in it, with the right cells, and spread over it).
- Single zone unchanged: seeded `gradb` and `gradb_gpu` are bit for bit
  the same as G1, `z10_gpu` has the same counts and velocity sums within
  2e-11 m/s. CPU and GPU `ctest` pass (except `MPITest`).

**G3. Multi-zone boundaries.** Done.

- `include/boundary_hd.h` (`FLAN_HD`) holds the per-particle boundary
  routines, taken from the loop body of `check_bounds_cpu`: `time_bc`,
  `check_faces` (absorbing, wall, core, `lcfs` and periodic faces,
  crossings into another zone up to `max_crossings`, the turn of a
  particle that goes around a toroidal wedge), `in_wall_cell` (the `chi`
  wall mask) and `turn_velocity`, with the helpers `absorbing_bc`,
  `periodic_bc`, `core_bc`, `lcfs_bc` and `core_face`. `core_face` uses the
  slot's RNG and picks among the face's core zones in the same order as
  before. The settings are a `BoundaryParams` (`Boundary::make_params`),
  the same for both. `check_bounds_cpu` calls them per particle and does
  the deposition and self-sputtering itself; `check_bounds_kernel` calls
  them per thread. The GPU's own copies of the helpers are gone.
- `ZonesView` has the core zones of each face (`core_offset`,
  `core_zones`), and `BackgroundDevice` has `chi` (when the background has
  one).
- What happened to the particles is counted per thread in a `FateTally`
  and added to `unsigned long long` counters on each GPU
  (`IonizRecombWarningsDevice::fates`) with `atomicAdd`, only where
  nonzero. `add_fates_from_device` adds them to the host's
  `Boundary::Fates`, and `follow_impurities` prints and returns them as on
  the CPU (`self_sputtered` stays 0 until G7). `Fates` moved to
  `boundary_hd.h` so the CUDA files can use it without MPI.
- The multi-zone, wedge and `lcfs` stops are gone.
- Tests: `test_boundary_zones_cuda` puts particles over every zone and up
  to 40% of its size past its faces (some past the last time) on the
  SOLEDGE3X test mesh (crossings, absorbing, wall and core faces, `chi`),
  the HERMES test grid (six zones, core faces, a wedge) and the slab with
  a core x face, periodic y and `lcfs` z. It runs `check_bounds_cpu` one
  particle at a time with slot i's GPU generator, `pcg32(seed + i, i)`, so
  state, zone, position and velocity must agree to 1e-12 even after a
  random core pick, and every fate count must be equal. With the GPU's
  velocity turn left out, the HERMES case fails. `test_boundary_cuda` now
  calls the shared helpers. The crossing cases of `test_cross_zones` are
  covered by these, since crossings are now the same code on both.
- Single zone unchanged: seeded `gradb` and `gradb_gpu` are bit for bit
  the same as G2, `z10_gpu` within 1e-14 (atomic sums). CPU and GPU
  `ctest` pass (except `MPITest`).
- The CPU's walls through the shared routines: `reg_self_sputter` (y05,
  y125, w_te20, w_te50), `reg_prompt_redep` (p01, p1) and `reg_drifts`
  pass their checkers.
- End to end, on the login node (one A100, 32 threads, shared): 4000 W
  from a point in the main SOL of a HERMES background, and spread over
  the near SOL zone of a SOLEDGE3X background, 2e-5 s, no sputtering, no
  sheath, thermal force off. The impurity content of each zone and its
  mean charge agree between CPU and GPU within the Monte Carlo error, as
  do the fates, and `check_hermes_bkg.py` passes on the GPU output. The
  main loop takes 2.5 min on the CPU and 8.7 s on the GPU (HERMES), and
  4.9 min and 16 s (SOLEDGE3X).

**G4. Thermal force.** Done.

- `include/collisions_hd.h` (`FLAN_HD`) holds the skew: `SkewCounts`,
  `SkewShape`, `skew_eps`, `skew_at` (eps and the unit vector along B
  from q_par and B, 0 where there's no heat flux or field),
  `skew_factor`, `skew_bound` and the sampler `sample_bkg` (the rejection
  loop with the |eps| cap). `skew_shape` and the push constants moved
  there too, host only. `nanbu_coll` and `nanbu_coll_kernel` both call
  them; the GPU's own unskewed sampler is gone, and draws the same
  random numbers as before.
- `BackgroundDevice` has `q_par_i` and `q_par_e` (when the background
  has them), interpolated along the field line like ne, with B from the
  stencil, as on the CPU. `OptionsDevice` has `imp_thermal_force`.
- Each GPU counts the sampler's draws per species
  (`IonizRecombWarningsDevice::skew`, `atomicAdd` where nonzero);
  `skew_counts_from_device` reads them and `add_skew_counts` adds them to
  the rank's, so `flan.cpp` reduces and prints them as before.
- The GPU stop in `check_thermal_force` is gone.
- Tests: `test_collisions_cuda` runs `collision_gpu` and `nanbu_coll` on
  the same particles (the CPU one at a time with slot i's GPU
  generator) on the slab with strong heat fluxes (eps 0.3, electrons
  sub-cycled) and on the HERMES test grid (conduction heat fluxes, six
  zones), for ions and electrons: velocities agree to 1e-9 of the speed,
  and every sampler count is equal. With the GPU skew turned off, both
  fail. `ThermalForce.Checks` now expects no stop on the GPU.
- `reg_thermal_force z10_gpu` passes `check_collisions.py` like `z10`
  (drift and coefficients within the statistical error of the CPU run,
  and the same share of draws clipped); its main loop takes 48 s
  against 10 min on 32 CPU threads.
- Seeded `gradb_gpu` is bit for bit the same as G3 and `z10_gpu`
  (friction) within 2e-14. CPU and GPU `ctest` pass (except `MPITest`).
- End to end, the G3 HERMES case with the thermal force on: zone
  contents, mean charges, fates and sampler counts agree between CPU
  and GPU within the Monte Carlo error, and `check_hermes_bkg.py` passes.
  The rejection sampling makes the GPU main loop about twice as long as
  without the thermal force (16 s against 8.7 s; 2.6 min on the CPU),
  since a warp waits for its slowest thread and the strongly skewed ions
  need about 7 draws per sample there.

**G5. Sputtered source and deposition.** Done.

- `include/sputtering_hd.h` (`FLAN_HD`) holds flat views of the surface
  elements, the yield table, the thermal quadrature, the source and the
  deposition tallies (`ElementsView`, `YieldView`, `Quadrature`,
  `SourceView`, `DepositionView`), and the routines `get_yield`,
  `sample_impact_energy`, `sample_thompson`, `cosine_direction`,
  `sample_birth`, `add_impact`, `add_lifetime_and_charge`,
  `add_into_wall_cell` and `add_at_wall_face`. `Surface`, `YieldTable`,
  `Source` and `Deposition` give views of their own data (`view()`) and
  call these, so the CPU and the GPU start atoms and tally their landings
  the same way. Where the CPU built a cumulative array to pick from
  (impact samples, quadrature nodes), the shared code forms the running
  sum twice instead, which picks the same entry and allocates nothing.
- `Deposition` keeps its wall cells and wall faces as sorted keys with
  their elements (searched by bisection) in place of hash maps, so the
  device searches the same arrays; the scalar tallies are one array. The
  tallies are added to through an adder passed to the shared routines:
  `omp atomic` on the CPU, `atomicAdd` on the GPU.
- `Source::to_device` copies the source's arrays to a GPU
  (`SourceDevice`, with the weight `source_rate / imp_num`);
  `fill_slots_gpu` starts atoms from it with `sample_birth`, neutral, with
  no prompt window. `Deposition::to_device` makes empty tallies on a GPU
  with the lookups and bins (`DepositionDevice`); `check_bounds_kernel`
  adds to them at wall faces and in wall cells, and
  `Deposition::add_from_device` adds them to the host's after the run,
  before the ranks' tallies are reduced. Each GPU has its own copies.
- `SlotsDevice` has `prompt_by` and `generation`. The ionization kernel
  takes `imp_iz_max_charge` and starts the prompt window on the first
  ionization from neutral, with |B| from the stencil, as on the CPU.
- The stops on the sputtered source, `imp_iz_max_charge` and `slab_wall`
  are gone; `sput_self` now stops on the GPU (G7).
- Tests: `test_sputtering_cuda` compares, with the CPU run one particle
  at a time on slot i's GPU generator: births from the HERMES test
  grid's targets (impact energy from the formula) and from the SOLEDGE3X
  test mesh's staircase (impact samples from the sheath), the ionization
  step with a charge limit and prompt windows, and the deposition tallies
  at the HERMES targets and in the SOLEDGE3X wall cells. Births and
  charges agree to rounding, the tallies to 1e-12. With any of the GPU
  routines turned off (births, charge limit, prompt window, either
  deposition routine), the tests covering it fail. The CPU results are unchanged
  (`test_sputtering_cpu` and the CPU wall regressions pass).
- `reg_prompt_redep` has `<variant>_gpu` for the cases without the
  sheath; `p01_gpu` and `p1_gpu` pass `check_redeposition.py` like the
  CPU runs (prompt fraction 0.978 and 0.590 on the GPU, 0.976 and 0.586
  on the CPU, against Fussmann's 0.977 and 0.587), with a main loop of
  about 1 s against 23 s and 86 s on 16 CPU threads. A `_gpu` variant of
  `reg_sheath_impacts` was left out: that case checks only the source,
  which is computed on the host for either path, and the births from
  impact samples are compared directly in `test_sputtering_cuda`.
- End to end, the HERMES targets as the source, sheath off, CPU against
  GPU: the source rate is identical, and the redeposited and prompt
  fractions, each zone's share of the redeposition, the fates and the
  mean charge, energy and angle on landing agree within the Monte Carlo
  error. Main loop 7 s on the GPU against 1.6 min on 32 CPU threads.
  The same with the SOLEDGE3X wall (staircase elements, atoms absorbed
  in wall cells): identical source, 98.8% and 98.7% redeposited, each
  zone's share within the Monte Carlo error; main loop 10.5 s on the
  GPU against 39 s on the CPU, where most atoms land within a
  microsecond. With nearly every atom coming back, the net erosion is a
  small difference of large numbers, and agrees only as well as the
  redeposition does (to about 0.005 of the source per zone).

**G6. Sheath field and sub-steps.** Done.

- `include/sheath_hd.h` (`FLAN_HD`) holds `Model`, `Params`, the
  profiles (`potential`, `field`, `range`), the sub-step count, and a
  `SheathView`: the wall elements as flat arrays (zone, face, the cell's
  extent across the face, thickness, normal), `Params` per frame and
  element, and the plasma cells next to the wall as sorted keys with
  each one's elements in CSR form, searched by bisection. `distance`,
  `field_at` and `substeps` on a view replace the bodies of
  `Sheath::Sheath`'s members, which now call them, so the CPU and the GPU
  find the same element and field. `Sheath::to_device` copies the arrays
  to a GPU (`SheathDevice`, freed with `free_sheath`); each GPU gets its
  own copy in `follow_impurities`.
- `cuda/boris.cu`: the push of one particle is a device function
  (`push_one`, the GPU's `push_one_cpu`) that adds the sheath field to
  the interpolated E. With the sheath on, `update_velocity_gpu` runs
  `count_substeps_kernel` (into the new `SlotsDevice::nsub`), the plain
  push of the rest, and `substep_kernel`, `Boris::substep_cpu` per
  thread: the half-kick, the sub-steps with the cell updated in the zone,
  and the kick to the crossing point for a particle that leaves the grid.
  `step_kernel` advances only the time of a sub-stepped particle. The
  half step back at the start adds the field without sub-steps, as on
  the CPU. `flan.cpp` builds the sheath for GPU runs too; the warning is
  gone.
- Decision 4's per-thread loop was kept. On the HERMES targets every
  sub-stepped particle took 2 sub-steps (about 150 to 230 particles per
  step of 4000), and the count, push and sub-step kernels together took
  about 8 s of a 20000-step run, measured with CUDA events on a GPU
  another process was using. No compaction pass is needed for this; G8
  profiles the full cases.
- Tests: `test_sheath_cuda` pushes 3000 particles (some neutral) placed
  in the cells next to the wall, on the slab wall with `"brooks"` and on
  the HERMES test grid's targets with `"borodkina"`, once with the field
  alone and once through the full near-wall step. Sub-step counts agree
  exactly, velocities and positions to 1e-8 after up to hundreds of
  sub-steps. With the GPU's sheath field or its sub-steps turned off,
  both tests fail. `test_sheath_cpu` and the CPU sheath regressions
  are unchanged.
- `reg_prompt_redep`'s `_gpu` suffix now takes the sheath cases. On the
  GPU, `drop0_gpu` passes `check_redeposition.py` (0.591, against 0.587
  on the CPU and Fussmann's 0.587); `brooks_a2`, `brooks_a5`,
  `borodkina_a2` and `borodkina_a5` agree with the CPU within the Monte
  Carlo error in prompt fraction, charge, impact energy and angle; so
  does `brooks_ref` (impact energy 104 eV, angle 13 degrees, on both).
  Each run took 2 to 23 s on the GPU against 20 to 127 s on 16 CPU
  threads. A `_gpu` variant of `reg_sheath_angles` was left out: it calls
  the push routines itself rather than running `flan()`, so it would need
  CUDA in the regression build; `test_sheath_cuda` compares the same
  kernels with the CPU on the same slab wall, and the `reg_prompt_redep`
  cases run the GPU's loop with the sheath.
- End to end, the HERMES targets as the source with `"borodkina"`, CPU
  against GPU: identical source rate, redeposited 0.682 and 0.672,
  prompt 0.678 and 0.670, each zone's share, the fates and the mean
  charge, energy and angle on landing all within the Monte Carlo error.
  On a compute node (one A100 against 64 CPU threads) the main loop took
  9.4 s on the GPU with the sheath, 6.5 s without it, and 31 s on the
  CPU with it. The sheath regressions `brooks_a5_gpu` and
  `borodkina_a5_gpu` took 3.2 s and 1.4 s there.

**G7. Self-sputtering.**

- `SelfSourceDevice`: the W-on-W yield table, the per-element parameters
  it needs from the source, the minimum weight, the maximum generation,
  and the tallies (`m_self`, `m_gen_weight`, `m_gen_landed`, the sums and
  counts) as device arrays and counters.
- `rebirth` as a shared function called from the boundary routine,
  with Russian roulette from the slot's RNG. A device flag for the
  avalanche (`sput_self_max_gen`), read after the step, stops the run as
  on the CPU.
- Regression: `_gpu` variants of `reg_self_sputter` (`y05`, `y125`).
  End to end: the HERMES target case with self-sputtering, CPU against
  GPU.

**G8. Performance and memory.**

- Profile a full HERMES and a SOLEDGE3X case with everything on (`nsys`,
  `ncu` on the push, collision and boundary kernels): occupancy, register
  use, divergence near the wall.
- Memory: the background, its frames and the stencils for large BOUT++
  runs on a 40 GB A100. If a background doesn't fit, stop with a message
  saying how much it needs (decision 5).
- Several GPUs per rank and several ranks: `MPITest` on a compute node,
  and a 1-node, 4-GPU run against 1 GPU.

**G9. Docs.**

- Remove "CPU only" from `input-options.md`, `running/soledge3x.md`,
  `running/hermes.md`, `physics/sputtering.md`, and the code comments.
- `code-structure.md`: the new device structs in the host/device table,
  the `FLAN_HD` shared functions, and how to add a feature that runs on
  both.
- `running-flan.md` GPU section: what runs on the GPU, memory limits, a
  note on the speed-up to expect.
- Update the "Later" notes of the HERMES and sheath plans to point here.

## Decisions

All seven were accepted as recommended (2026-11-24).

1. **One path or a separate single-zone fast path.** Recommended: one
   path (G1), with constant memory removed, kept only if G1's timing
   shows single-zone runs get noticeably slower.
2. **Shared `FLAN_HD` code or separate CUDA copies.** Recommended:
   shared functions for everything new and for what each phase touches,
   leaving untouched duplicates (such as `nanbu_s_a.cuh`) for later.
3. **Where births are sampled.** Recommended: on the device, from a
   device copy of the source (G5). Sampling on the host would need a
   copy to and from the GPU every step.
4. **Near-wall sub-steps.** Recommended: per-thread loop first, a
   compacted near-wall kernel only if profiling shows the warps waiting
   on it.
5. **Backgrounds too big for one GPU.** Recommended: stop with a clear
   message for now. Keeping only the frames around the particles' times
   on the device is a larger change, left for later.
6. **Order.** Recommended: G0 to G3 first, since they make SOLEDGE3X
   and HERMES runs work on the GPU at all and everything later depends
   on the zones; then G4, G5, G6, G7 in that order (each depends on the
   one before, except G4).
7. **How close CPU and GPU must agree.** Recommended: deterministic
   tests to rounding; regression and end-to-end results within the
   Monte Carlo error (e.g. 3 sigma), not bit for bit.
