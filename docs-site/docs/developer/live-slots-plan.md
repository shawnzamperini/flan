# Live-Slot List Plan

This page is the development plan for running the GPU kernels over a list
of the slots that hold a live particle, instead of over every slot. It is
the first item of "Later" in the [GPU Plan](gpu-plan.md#later). The work
is on the branch `gpu-live-slots`.

## Why

Every GPU kernel in the main loop runs one thread per slot, and a thread
whose slot is dead returns at once. When a run has fewer particles than
`slot_cap` (about a million by default), every particle gets a slot at
the start and slots aren't refilled, so as particles die the live ones end
up scattered among dead slots. A warp then runs with only a few of its 32
lanes doing work, but takes as long as its slowest lane.

G8 measured this on a HERMES case with everything on: about 30% of the
slots live mid-run and 3% near the end, and in the collision kernel (61%
of the GPU time) only 2 to 3 of 32 lanes busy per instruction. The
collision sub-cycles weren't the cause (almost every particle has one).
Packing the live slots into the first threads should give about 1.5 to 2
times in the busy part of a run. It won't help the time a step takes when
only a few particles are left, which is the latency of one particle
through the kernels.

## Approach

### The list

`fill_slots_gpu` already lists the dead slots each step
(`cub::DeviceSelect::If` into `SlotsDevice::dead_indices`) and copies
their number to the host, which is the step's one wait for the GPU before
"Checking particle state". After the fill kernel, a second selection
lists the live slots in increasing order into a new
`SlotsDevice::live_indices`. The list is in slot order, the same every
run, so seeded runs stay reproducible.

The host already knows how long the list is: `alive_slots`, the number of
live slots after the fill, is worked out from the dead count. Nothing
between the fill and the next step's `step_kernel` kills a particle (only
`step_kernel`, for `imp_max_time`, and `check_bounds_kernel` do), so
`alive_slots` stays exact until the list is rebuilt. The launches are
sized from it, with no new copy from the device.

Particles that die during a step (in `step_kernel` or the boundary check)
stay in the list until the next fill. Each kernel keeps its check of the
slot's state, as now, so they're skipped.

### The kernels

Each per-particle kernel takes its slot from the list:

```cpp
const int k {static_cast<int>(blockIdx.x * blockDim.x + threadIdx.x)};
if (k >= slots_d.num_live) return;
const int i {slots_d.live_indices != nullptr ? slots_d.live_indices[k] : k};
```

behind one small device helper (`Slots::slot_of_thread`), and launches
`ceil(num_live / 256)` blocks. Those kernels are: ionization and
recombination, the two collision launches, the Boris push, the sub-step
count and the sub-steps, the step, the boundary check, finding the cell
and recording the statistics. `fill_slots_kernel` already runs over the
dead list and `init_rngs` over every slot; they don't change.

Each slot keeps its own random number generator (`rngs_d[i]`) and all
the per-particle results are written to slot i, so which thread handles
a slot doesn't change anything: the trajectories and counts are the same
bit for bit, and the sums made with `atomicAdd` (statistics, deposition)
differ only by rounding, as between two runs today.

### Without a list

The unit tests call the `*_gpu` functions on slots they set up
themselves, without a fill. With no list (`live_indices == nullptr`, the
default) and `num_live = N`, the kernels run over every slot as now, so
the tests work unchanged and a test can compare the two ways.

### The CPU

Unchanged. A dead slot costs a CPU thread only the check of its state,
and the thread moves straight on to the next slot; it doesn't wait for
the other slots of a warp the way a GPU lane does.

## Phases

**L0. Baseline.**

- Seeded reference outputs from `main` for the end-to-end cases (HERMES
  `all`, SOLEDGE3X `self`, a Gkeyll case) and the `_gpu` regression
  variants, to compare L1 and L2 against.
- Temporary instrumentation, not committed: how many of the live slots
  are neutral each step, which decides decision 1.

**L1. The live list.**

- `SlotsDevice::live_indices` (allocated and freed with the others) and
  `num_live`; `fill_slots_gpu` builds the list and sets `num_live`
  after filling; `Slots::slot_of_thread`; every kernel listed above uses
  it and is launched over `num_live` threads.
- Tests: a `*_cuda` test that the list holds exactly the live slots, in
  order, after fills with dead, live and refilled slots; and that a
  kernel run over the list gives the same slots as one run over every
  slot (push, collisions, boundary).
- Checks: `ctest` on both builds; the seeded end-to-end runs equal to
  L0's (counts and trajectories exactly, sums to rounding,
  `g0/cmp.py`); the `_gpu` regression variants pass.

**L2. The all-dead check from the host's count.**

- With `alive_slots` exact after the fill, the run is over when it is 0
  and no particles are left to start. That replaces `all_dead_kernel`
  and its copy from the device, one wait for the GPU each step (decision
  3). Same checks as L1.

**L3. A charged-only list for collisions** (only if L0 shows many live
slots are neutral, decision 1).

- Ionization runs just before the collisions and changes the charges,
  so the list is built on the device after it, and its length stays on
  the device: the collision kernels are launched over `num_live` threads
  and a thread past the device-side count returns.

**L4. Measure.**

- Compute-node job: the G8 cases against L0 on the same node, 1 and 4
  GPUs, nsys for the kernel times and ncu for the busy lanes. Try blocks
  of 128 threads (decision 4).
- Record the results here and in the GPU plan's "Later".

**L5. Docs and merge.**

- `code-structure.md`: how a kernel finds its slot (`slot_of_thread`),
  for new kernels. Then merge to `main`.

## Decisions

All five were accepted as recommended (2026-11-26).

1. **Neutral particles in the list.** Recommended: one list of every
   live slot, used by all the kernels (ionization and the push need the
   neutrals), and a second, charged-only list for the collisions (L3)
   only if L0 shows enough neutrals to make it worth a selection each
   step.
2. **Kernels without a list.** Recommended: fall back to every slot when
   there's no list, so the unit tests and anything calling the `*_gpu`
   functions directly keep working, rather than making every caller
   build one.
3. **The all-dead check.** Recommended: use the host's count (L2), since
   it is exact after the fill, and drop the kernel and the wait it costs
   each step. The alternative is to leave it as is and keep L1 the only
   change.
4. **Block size.** The kernels use 133 to 255 registers per thread, so
   an SM fits one block of 256 threads, and a launch that's a little
   over a multiple of the SM count leaves most SMs idle for its last
   blocks. Recommended: measure 128 in L4 and keep it only where it's
   faster; otherwise 256 as now.
5. **How close the results must be.** Recommended: the same as within
   one GPU run today: counts and trajectories exactly, `atomicAdd` sums
   to rounding.
