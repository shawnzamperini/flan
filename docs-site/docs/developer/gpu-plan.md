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
| More than one zone (SOLEDGE3X, HERMES) | Stops | `Background::to_device` calls `check_single_zone` |
| Toroidal wedge (velocity rotation across the z ends) | Stops | `Background::to_device` |
| Field-aligned (zShift) interpolation, HERMES | Not reachable (multi-zone) | `Utilities::interp_aligned`, `push_one_cpu` |
| Ghost cells, zone crossings, `face_core`, `face_wall`, the `chi` wall mask | Not reachable (multi-zone) | `Boundary::check_bounds_cpu` |
| `lcfs` z boundary (Gkeyll) | Stops | `crosscheck_options` |
| Thermal force (heat-flux-skewed sampling) | Stops if the background has heat fluxes | `check_thermal_force` |
| Sputtered source (`imp_source_opt = "sputter"`) | Stops | `flan.cpp` |
| Deposition and prompt tallies, impact histograms | Not reachable (needs the source) | `Sputtering::Deposition` |
| Sheath field and near-wall sub-steps | Left out, with a warning | `flan.cpp`, `Sheath::field_at`, `Boris::substep_cpu` |
| Self-sputtering | Not reachable (needs the source) | `Sputtering::SelfSource::rebirth` |
| `imp_iz_max_charge` (test option) | Stops | `flan.cpp` |
| Boundary fates (`Boundary::Fates`) | All zeros | `main_loop` |
| `Slots` fields `zone`, `prompt_by`, `generation` | Host only | `slots.h` |

Things that are host-side setup and need nothing on the GPU: reading the
backgrounds (Gkeyll, SOLEDGE3X, both HERMES formats), the field-aligned
derivatives, the conductive heat fluxes and the inductive E_par (all
computed into the background before the run), `Sputtering::Surface`, the
source's rates and the D+ impact table (`Sheath::ImpactTable`, built at
startup), and saving the results.

## What the GPU has now

- Device structs of raw pointers for the background, slots, statistics,
  OpenADAS and options (see [Code Structure](code-structure.md)). The
  background's 1-D grids (times, cell centers and edges) are in constant
  memory, `MAX_1D = 500` each, so only one background fits on a GPU.
- The background data on the device is a single zone with no ghost cells,
  which has the same layout as a `Vector4D`. A `ZonedVector` is already
  flat (zone after zone, with per-zone offsets and strides), so the
  multi-zone data can be copied as it is.
- Interpolation through an 8- or 16-point stencil built once per particle
  per step (`include/cuda/interpolate.cuh`), clamped at the grid's edges.
- Kernels for every main loop step, with the Nanbu collisions (including
  sub-cycling), OpenADAS ionization and recombination, the Boris push in
  curvilinear coordinates, absorbing/periodic/core boundaries, filling
  slots from the `imp_*start_opt` options, and the statistics.
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

**G0. Baseline.**

- Build `build_gpu` from the branch; all `*_cuda` tests and the existing
  `z10_gpu` and `gradb_gpu` variants pass.
- Time a single-zone GPU run (slab and Gkeyll) by main loop section
  (`Timer`), as the reference for G1.
- Add the `FLAN_HD` macro (`include/host_device.h`) and move one small
  routine over (e.g. `trilinear_interpolate`) to check the build in both
  directories and the docs' Doxygen preprocessing (`PREDEFINED` in
  `mkdocs.yml`).

**G1. Zones on the device.**

- `Zones::ZonesDevice` (`include/zones_device.h`): per-zone nx, ny, nz,
  the ghost width, the offsets of each zone's centers and edges in flat
  arrays, the face codes (`num_faces` per zone), the toroidal, Cartesian
  and wedge flags, zShift offsets and table, and the `ZonedVector`
  offsets and strides (all fields share one layout per background, so
  one set). `Zones::to_device`/`free_zones`.
- `BackgroundDevice` holds a `ZonesDevice` and the flat `ZonedVector`
  data, ghosts included. Remove `check_single_zone` and the constant
  memory from `Background::to_device`.
- `SlotsDevice` gains `zone`.
- Shared `FLAN_HD` versions of `Zones::center`, `neighbor_index` and
  `bracket` (index finding with ghosts), used by the CPU too.
- Every kernel takes its zone from the slot. Single-zone behaviour is
  unchanged: GPU runs of the slab and Gkeyll cases give the same results
  as on main (bit for bit if the arithmetic didn't change order, otherwise
  within rounding), and run no slower than G0.
- Tests: `test_zones_cuda` (centers, neighbor indices and brackets on the
  SOLEDGE3X and HERMES test grids against the CPU), `test_interpolate_cuda`
  extended to several zones with ghosts.

**G2. Multi-zone push, collisions, ionization and statistics.**

- Interpolation stencils per zone, using ghost cells at zone edges as
  the CPU does.
- Field-aligned interpolation on the device: a shared version of
  `Utilities::interp_aligned` (scalars, and vectors rotated per cell
  frame), used in the Boris push where the zone has a zShift.
- Boris (`interp_recp`, `interp_tang` per zone), collisions and
  ionization/recombination with the zone.
- `StatisticsDevice` per zone (flat `ZonedVector` layout), and
  `fill_slots_gpu` starting in `imp_zone_start_val`.
- Tests: the push of one particle on a HERMES test grid against
  `push_one_cpu`, aligned and not (`test_boris_cuda`).

**G3. Multi-zone boundaries.**

- A shared per-particle boundary routine (from the loop body of
  `check_bounds_cpu`): absorbing and wall faces, crossings into another
  zone (`max_crossings`), periodic faces, the wedge velocity rotation,
  `core_face` (with the slot's RNG), the `lcfs` z boundary, and the `chi`
  wall mask. The CPU loop calls it per particle; `check_bounds_kernel`
  calls it per thread. Deposition hooks are left as no-ops until G5.
- `Boundary::Fates` counted on the device (one counter per field,
  `atomicAdd` on `unsigned long long`), summed into the host fates.
- Remove the wedge and `lcfs` guards.
- Tests: `test_boundary_cuda` with a list of positions crossing every
  kind of face on the SOLEDGE3X and HERMES test grids and on a wedge,
  against the CPU routine; `test_cross_zones` cases on the device.
- End to end: `s3x_example` and `hermes_example` (point source, no
  sputtering, thermal force off) on CPU and GPU; density by zone within
  the Monte Carlo error, fates agree, `check_hermes_bkg.py` passes on the
  GPU output. This is the first point where a SOLEDGE3X or HERMES run
  works on the GPU.

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
