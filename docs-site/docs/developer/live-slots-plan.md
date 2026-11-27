# Live-Slot List Plan

This page is the development plan for running the GPU kernels over a list
of the slots that hold a live particle, instead of over every slot. It is
the first item of "Later" in the [GPU Plan](gpu-plan.md#later). The work
was done on the branch `gpu-live-slots` and merged to `main` on
2026-11-27.

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

Done (2026-11-26). The references are in `$SCRATCH/flan_gpu/live/l0`
(`refs.sh` runs them): end-to-end HERMES `all` and `sput`, SOLEDGE3X
`self`, the Gkeyll timing case and the seeded friction case. Neutral
share of the live slots, sampled every 1000 steps after the first:

| Case | Neutral share of live slots |
|---|---|
| HERMES, sputtered source | 46% at step 5000, 26% at 10000, 16% at 15000, 10% at 20000; 37% of the live particle steps |
| SOLEDGE3X, sputtered source | All neutral at birth, but nearly all gone by step 2000; few neutrals after |
| Gkeyll, started charged | None |

So in a sputtered-source run, a third or more of the live slots early on
are neutral, falling as they ionize. Whether L3 is worth doing is left to
L4's measurement of the busy lanes in the collision kernel.

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

Done. The selection reuses the dead list's work space; `live_indices`
and its device-side count are allocated by the first fill. Every
per-slot kernel in `src/cuda/` takes its slot from `slot_of_thread` and
is launched over `num_blocks` blocks (none when nothing is live).
`LiveSlotsCUDATest` checks that there's no list before a fill, that the
list is the live slots in order after fills that kill and refill
scattered slots, and that collisions, the push, the step and the cell
search give bit-identical slots with and without it. Two mutations
(ignoring the list, listing the dead slots) each fail it. Against L0,
the seeded runs differ only in `atomicAdd` sums (the mean velocities,
charge, density and redeposited energy), by at most 1e-15 relative;
every count, and the printed summaries, are the same. The `_gpu`
regression variants pass their checkers (collisions, drifts,
redeposition, self-sputtering); `ctest` passes on both builds except
`MPITest` (compute node only).

**L2. The all-dead check from the host's count.**

- With `alive_slots` exact after the fill, the run is over when it is 0
  and no particles are left to start. That replaces `all_dead_kernel`
  and its copy from the device, one wait for the GPU each step (decision
  3). Same checks as L1.

Done, with L1: `all_dead_wrapper` returns `num_live == 0` once there's
a list, and falls back to `all_dead_gpu` before the first fill. Checked
with L1's runs.

**L3. A charged-only list for collisions** (only if L0 shows many live
slots are neutral, decision 1).

- Ionization runs just before the collisions and changes the charges,
  so the list is built on the device after it, and its length stays on
  the device: the collision kernels are launched over `num_live` threads
  and a thread past the device-side count returns.

Dropped after L4: the busy lanes in the collision kernel are set by
divergence inside it, not by which slots share a warp (see L4), so
removing the neutrals from its list can't gain much.

**L4. Measure.**

- Compute-node job: the G8 cases against L0 on the same node, 1 and 4
  GPUs, nsys for the kernel times and ncu for the busy lanes. Try blocks
  of 128 threads (decision 4).
- Record the results here and in the GPU plan's "Later".

Done (2026-11-26), jobs in `$SCRATCH/flan_gpu/live/l4`. `main`, this
branch and this branch with 128-thread blocks were run one after another
on one node. Wall time (s), HERMES `all` unless noted:

| Case | `main` | list | list, 128 |
|---|---|---|---|
| 1 GPU, 40000 | 84.3 | 78.0 | 74.6 |
| 1 GPU, 160000 | 155.4 | 113.3 | 108.0 |
| 4 GPUs, 160000 | 83.3 | 76.5 | |
| Gkeyll, 1 GPU | 9.2 | 8.7 | 7.6 |

With few particles the list is slower. Repeats (three of SOLEDGE3X
`self`, two of HERMES at 4000, alternating) put "Following impurities"
at 10.25 s for `main` against 11.27 s with the list for SOLEDGE3X, and
28.2 s against 29.5 s for HERMES. The wall times of these runs vary more
than that, from reading the background. nsys on SOLEDGE3X: the collision
kernel takes 20% longer with the list and the sub-steps 10%; the others
are the same.

ncu at step 8000 of HERMES at 40000, where 73% of the slots are live
(114 blocks instead of 157): busy lanes per instruction in the collision
kernel go from 2.3 and 3.0 (ions, electrons) to 2.7 and 3.8, the
sub-steps from 4.0 to 5.3, the push from 15.6 to 21.4 and ionization
from 17.2 to 23.5. Across the run nsys has the collision kernel 8%
faster.

So the idle lanes in the collision kernel come from divergence inside
it (the rejection loop of the skewed sampler, sub-cycles, neutrals
returning early), not from dead slots. A warp runs its lanes' different
paths one after another, so packing particles into fewer warps only
helps when there are more warps than the GPU keeps resident. With fewer
particles than that, `main` gives each its own warp, on its own SM, and
they run in parallel; the list puts them in a few warps that run them in
turn.

128-thread blocks (decision 4), ncu: ionization 85 to 51 us, collisions
about 5% faster, the push 2%, the sub-steps the same. Kept for those
(`Slots::heavy_block_size`); the other kernels stay at 256. Results are
unchanged apart from the rounding of `atomicAdd` sums.

**L4b. Spreading the list over the warps** (added after L4).

- While there are fewer live slots than the GPU has resident warps for
  the heavy kernels, give each warp only as many as fills the GPU, down
  to one (`Slots::per_warp_for`, `Slots::list_position`); pack them 32
  to a warp once there are more. `resident_warps` is set by the first
  fill from the SM count and the registers per SM, taking the heavy
  kernels at 255 registers per thread.

Done (2026-11-27). `LiveSlotsCPUTest` checks that each live slot goes to
exactly one thread with no more than `per_warp` to a warp (dropping the
lane cut-off fails it), and `LiveSlotsCUDATest` that its comparison case
is spread. The seeded runs match L2 to rounding and the logs are the
same. On one node, alternating, "Following impurities" (s):

| Case | `main` | packed (L4) | spread |
|---|---|---|---|
| SOLEDGE3X `self` (2 runs each) | 10.26 | 11.19 | 9.82 |
| HERMES, 1 GPU, 4000 | 28.2 | 29.3 | 25.0 |
| HERMES, 1 GPU, 40000 | 51.1 | 43.7 | 40.3 |
| HERMES, 1 GPU, 160000 | 125.4 | 84.6 | 80.4 |
| HERMES, 4 GPUs, 160000 | 56.5 | | 45.3 |
| Gkeyll, 1 GPU | 2.87 | 2.38 | 2.36 |

So the spread list is faster than `main` in every case, by 4% when only
a few particles are left for most of the run and by a third in the
largest. The rest of the collision kernel's idle lanes are divergence
inside it, which is left for later (GPU plan).

**L5. Docs and merge.**

- `code-structure.md`: how a kernel finds its slot (`slot_of_thread`),
  for new kernels. Then merge to `main`.

Done: `code-structure.md` covers `slot_of_thread`, the spreading and
`heavy_block_size`, and the GPU plan's "Later" points here.

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
