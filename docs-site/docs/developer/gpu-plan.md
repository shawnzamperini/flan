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
| Thermal force (heat-flux-skewed sampling) | Stops if the background has heat fluxes | `check_thermal_force` |
| Sputtered source (`imp_source_opt = "sputter"`) | Stops | `flan.cpp` |
| Deposition and prompt tallies, impact histograms | Not reachable (needs the source) | `Sputtering::Deposition` |
| Sheath field and near-wall sub-steps | Left out, with a warning | `flan.cpp`, `Sheath::field_at`, `Boris::substep_cpu` |
| Self-sputtering | Not reachable (needs the source) | `Sputtering::SelfSource::rebirth` |
| `imp_iz_max_charge` (test option) | Stops | `flan.cpp` |
| Boundary fates: `self_sputtered` | 0 (no self-sputtering) | `follow_impurities` |
| `Slots` fields `prompt_by`, `generation` (`zone` is on the device since G1) | Host only | `slots.h` |

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
- One host thread per GPU, slots split between ranks and then GPUs.
- Tests `*_cuda`, and the regression variants `z10_gpu` (friction) and
  `gradb_gpu` (drifts).

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

**G4. Thermal force.**

- The skewed sampler (`skew_factor`, `skew_bound`, the rejection loop
  with the |eps| cap) as shared functions, used by `collision_cpu` and
  the collision kernel. The heat fluxes go to the device with the
  background. The skew counts become device counters, reduced like the
  CPU ones.
- Remove the GPU stop in `check_thermal_force`.
- Regression: `reg_thermal_force` gets a `_gpu` variant; the HERMES
  example with the thermal force on, CPU against GPU.

**G5. Sputtered source and deposition.**

- `Sputtering::SourceDevice`: the CDF over (frame, element), each
  element's zone, cell bounds, normal, impact angle and the frame
  windows, the yield table's Thompson parameters, and for
  `sput_impact_opt = "sheath"` the impact table's samples flattened (one
  array of samples, offsets per point). `sample_birth` becomes a shared
  function, so the CPU and GPU draw births the same way. `fill_slots_gpu`
  uses it when a source is given, with the weight `source_rate / imp_num`.
- `SlotsDevice` gains `prompt_by` and `generation`.
- `Sputtering::DepositionDevice`: per-element rates, prompt rates and
  charges, impact energy and angle histograms, added to with
  `atomicAdd` in the boundary routine (`add_at_wall_face`,
  `add_into_wall_cell` as shared functions writing into plain arrays).
  Copied back and added to the host `Deposition` after the run.
- `imp_iz_max_charge` on the device (one option field and a check in the
  ionization kernel), so the prompt redeposition cases run on the GPU.
- Remove the stops in `flan.cpp`.
- Regressions: `_gpu` variants of `reg_prompt_redep` (sheath off) and
  `reg_sheath_impacts`. End to end: the HERMES target source, sheath off,
  CPU against GPU (source rate identical, redeposition and net erosion
  within the Monte Carlo error).

**G6. Sheath field and sub-steps.**

- `Sheath::SheathDevice`: the elements (normal, face point, thickness),
  `Params` per frame and element, and the cell-to-elements map flattened
  from the `unordered_map` to CSR (a sorted key array searched by
  binary search, or a per-cell offset array over the zones' cells if
  memory allows; G0's sizes decide).
- `potential`, `field`, `Sheath::distance`, `field_at` and `substeps` as
  shared functions. The Boris kernel adds the field; the step kernel and
  a device version of `Boris::substep_cpu` run the sub-steps per thread.
- A particle near the wall can take up to `max_substeps` (20000)
  sub-steps, and its warp waits for it. First do the straightforward
  per-thread loop and measure (G0's timer plus `nsys`). If the near-wall
  particles dominate, move them to a separate kernel after a compaction
  pass (decision 4).
- Remove the warning in `flan.cpp`; the GPU runs with `sheath_model` as
  set.
- Regressions: `_gpu` variants of `reg_sheath_angles` and
  `reg_prompt_redep` with the sheath. End to end: the HERMES target case
  with the sheath, CPU against GPU, and the run time.

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
