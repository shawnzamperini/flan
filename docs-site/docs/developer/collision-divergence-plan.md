# Collision Kernel Divergence Plan

This page is the development plan for keeping more lanes busy in the GPU
collision kernel (`nanbu_coll_kernel`, `src/cuda/collisions.cu`). It is
the "Divergence in the collision kernel" item of "Later" in the
[GPU Plan](gpu-plan.md#later). The work goes on the branch
`gpu-coll-divergence`, from `main`.

## Why

The collision kernel is the largest share of the GPU's time (61% in
G8's HERMES case, with everything on). G8 and the live-slot plan's L4
measured with ncu that only 2.3 to 3.8 of a warp's 32 lanes are busy per
instruction in it. The live list packs the warps with live slots, yet
the busy lanes barely changed (2.3 to 2.7 for ions, 3.0 to 3.8 for
electrons, with 73% of the slots live). So the idle lanes come from the
lanes of a warp taking different paths inside the kernel, which the warp
runs one after another. The Boris push and ionization, which branch
less, run 16 to 23 lanes.

Up to now the causes have only been guessed at: the rejection loop of the
skewed sampler, the collision sub-cycles and neutrals returning early.
Which of them matters is the first thing to find out, because the fixes
differ.

### What each cause could cost

These estimates were written before C0, which measured them (see C0
under Phases). Two were wrong: the sub-cycles cost nothing in HERMES,
and the sampler's bound M is not "at most about 2" but up to 33 for
the ions and 60 for the electrons, so it is the largest cost.

ncu's busy lanes per instruction is weighted by instructions, so a few
long paths count for more than many short ones. A rough count for each
cause:

- **Sub-cycles.** A particle gets `n_sub` sub-steps, from 1 up to
  `imp_coll_subcycle_cap` (128 by default), and a warp loops as many
  times as its lane with the most. G8 found that almost every particle
  has 1 and a few have 2. But if one lane in a warp has 100 while the
  others have 1, the warp does 100 sub-steps with one lane busy for 99 of
  them, and the warp averages about 1.3 busy lanes. A few such particles
  per warp would explain the measurement alone. Particles in cold, dense
  plasma near the targets are where to look.
- **The skewed sampler** (`sample_bkg`, `collisions_hd.h`). Each draw is
  accepted with probability F / M, and M is at most about 2, so a lane
  takes 1 to 2 draws on average. The warp waits for the lane with the
  most, which is a few draws at most. That costs a factor of 2 or 3 on
  the sampler, not 10 on the kernel.
- **The cached normal** (`pcg32::normal`). Every other call makes two
  values with Box-Muller (log, sqrt, sin, cos) and the next returns the
  cached one. Once lanes have taken different numbers of draws, their
  parity differs, and the warp runs both branches. That costs at most a
  factor of 2 on the draws.
- **Neutrals.** They return at once. They are 10% to 46% of the live
  slots in a sputtered-source HERMES run (L0), so they can cost up to
  about a factor of 2, early in the run.
- **Interpolation.** `interp_plasma` branches on whether the zone has a
  zShift table, so a warp with slots in two kinds of zone runs both.
  This only happens where zones meet.
- **The A(s) and chi branches.** `nanbu_calc_A_cuda` has three branches
  (one of them a table lookup). `nanbu_calc_chi_cuda` already works out
  every branch and then picks one.

## Approach

### Measure first (C0)

Build `build_gpu` as `RelWithDebInfo` (which already adds `-lineinfo`),
and profile the collision kernel at the same steps L4 used, with ncu's
source counters (`--section SourceCounters`, or `--set full` and the
source page). These give the warp and thread instructions for each
source line, and so where the idle lanes are. Also add temporary counters
(not committed) for each warp: the largest and total `n_sub`, the
largest and total sampler draws, and the number of neutral lanes. Their
ratios give each cause's cost directly.

The results decide which of C2 to C4 are done. Each is worth doing only
if C0 shows that what it removes is a large share of the idle lanes.

### Split the kernel (C2)

The ion and electron collisions are two launches of the same kernel, one
after the other. Each builds the interpolation stencil and interpolates
ne, Te, the flow, B and q_par again. A prepare kernel that does all of
this once per step for both species would store, for each slot:

- ne, Te, Ti, the flow (3), B (3), q_par_e and q_par_i (11 doubles), and
- `n_sub` for the ions and for the electrons (2 ints), 0 for neutral and
  dead slots.

The collision kernels then read these values and do only the sub-steps.
They work out `s_coeff`, eps and b from the stored values with the same
arithmetic as now. The interpolation then happens once instead of twice,
the collision kernels need fewer registers, and C3 has the `n_sub`
values it sorts by.

No random numbers are drawn in the prepare kernel, and each slot still
draws from its own generator in the same order (ions, then electrons).
So the results should stay the same bit for bit, as with the live list.
One thing could break that: the compiler may fuse multiplies and adds
differently once the code is split. If it does, the change is to
rounding only, and the seeded runs are compared as in L1, with sums to
rounding.

Memory: 96 bytes a slot, about 100 MB at the default `slot_cap`. The
arrays are allocated with the other slot arrays and covered by G8's
memory check.

### Order the work (C3)

Divergence from the sub-cycles goes away if lanes with similar `n_sub`
share a warp. After the prepare kernel, a stable sort of the live list
by a small key gives a collision list ordered by work. The key is the
bucket of `n_sub` (log2, 0 for neutral, so about 4 bits; ions and
electrons sorted separately, or by their sum, as decided in C0). This is
`cub::DeviceRadixSort::SortPairs` with `begin_bit`/`end_bit` limited to
the key. The sort is stable, so the order within a bucket is slot order,
the same every run. The neutrals come first (key 0) and are skipped:
the collision kernels start at the device-side count of neutrals and
read the list length on the device. That is L3's approach, and it adds
no wait for the host.

Which thread handles a slot doesn't change the slot's results, so this
step is bit-identical too. It costs one small sort a step, tens of
microseconds for 10^5 slots. When few slots are live and the list is
spread over the warps (L4b), the work order matters less, because each
warp has only a few lanes. In that case the sort is skipped (decision
3).

### Restructure the sampler (C4, only if C0 asks for it)

Written before C0. C0 found the size of M matters more than the
lanes waiting for each other, so C4 is now the tighter envelope below;
these two are kept for later, if C5 shows the waiting still matters.

If the sampler's draws are a large share after C3, make its loop the
same for every lane of a warp. Two candidates:

- **A fixed number of candidates per pass.** Draw K candidates at once
  and take the first one that is accepted, with another pass only if
  none is. With M at most 2, K = 2 or 3 makes the second pass rare.
- **Normals in pairs.** Make the three normals of a draw from two
  Box-Muller pairs, without the cached value. Every lane then runs the
  same transcendental functions, at the cost of throwing one value away.

Both change which random numbers are drawn, so the results change
statistically, not bit for bit. Accepted candidates have the same
distribution: taking the first accepted of K independent candidates is
still rejection sampling. `sample_bkg` is shared with the CPU, so the CPU
changes with it, as `CLAUDE.md` asks. The checks are the regression
collisions and thermal-force cases with their checkers (including the
`_gpu` variants), and the sampler's counts against `test_thermal_force`.

### A tighter envelope (C4, after C0)

C0 showed that the sampler's cost is the size of M, not only the
lanes waiting for each other: with M up to 33 or 60, even a warp of
lanes that all drew together would take tens of draws. Drawing K
candidates per pass doesn't change that. A tighter envelope does.

Since |w_par| <= w, the skew factor is bounded at every speed by

  F <= 1 + |eps| (|a| w + |b| w^3 + |c| w^5) = G(w).

f_M G is a mixture of four isotropic distributions: the Maxwellian
(weight 1), and w^k f_M for k = 1, 3, 5, with weights |eps| |a| <w>,
|eps| |b| <w^3> and |eps| |c| <w^5> (moments of the Maxwellian's speed:
<w> = 1.60, <w^3> = 6.38, <w^5> = 38.3). Each is easy to draw: an
isotropic direction and a speed with w^2 from a gamma distribution of
shape (k + 3) / 2 and scale 2. A candidate drawn from the mixture is
accepted with probability max(F, 0) / G(w). The mean number of draws
is the mixture's total weight, 1 + |eps| (|a| <w> + |b| <w^3> + |c|
<w^5>): at |eps| = 1, 3.0 for the ion shape and 3.1 for the electron
shape, against 33 and 60 now. It also removes the cut at 5 v_th, since
G bounds F at every speed: no draws above the bound, and the sample is
exact apart from the clipping of F < 0, which is unchanged.

The shapes for k = 1, 3, 5 are 2, 3 and 4, all integers, so each w^2
is a sum of 2, 3 or 4 exponentials, -2 ln U each. The direction is a
uniform cos(theta) and azimuth. The Maxwellian component is drawn as
now, from three normals.

Checked in Python (`envelope_check.py` and `envelope_check2.py` in the
C0 directory, 16 million candidates for each case, against the moments
of max(F, 0) f_M by quadrature). Mean draws per sample, current bound
and envelope:

| Shape | eps = 0.1 | eps = 0.5 | eps = 1 |
|---|---|---|---|
| Ions | 4.2 and 1.2 | 17.1 and 2.0 | 32.7 and 3.0 |
| Electrons | 6.9 and 1.2 | 30.2 and 2.0 | 57.9 and 3.0 |

The envelope's heat flux and <w_par^3> agree with the quadrature
within the statistical error in every case. The current sampler's are
3% low for the electron shape at eps = 1 (q / eps 0.493 against 0.508),
from the cut at 5 v_th, since the electron shape puts its lean on fast
particles.

The fixed-candidates idea above can still follow if C5 shows the lanes
waiting for each other matters once the draws are few.

### Sub-cycles that are still long

If C0 shows a few particles with `n_sub` near the cap, sorting puts them
together in their own warps, but each warp still runs them one after
another. That is the latency of one particle, and no reordering removes
it. Making those sub-steps cheaper (the cap, the target, or an
analytic limit for s >> 1) changes the physics and belongs in its own
plan. C0 records how often it happens.

## Phases

**C0. Measure.**

- Compute-node job, with the G8/L4 scripts (`$SCRATCH/flan_gpu/live/l4`):
  HERMES `all` at 40000 and 160000 particles, SOLEDGE3X `self`, Gkeyll.
  ncu source counters for the collision kernel at step 8000 (as in L4)
  and at a late step, and the temporary per-warp counters over the run.
- Record here how the idle lanes split between sub-cycles, the sampler,
  the cached normal, neutrals, interpolation and the rest, and the
  distribution of `n_sub` (and its largest value) for ions and
  electrons. Decide which of C2 to C4 go ahead.

Done (2026-11-27), job in `$SCRATCH/flan_gpu/coll/c0` (`c0.sbatch`; the
counters are a patch to a copy of the source in `src/`, and
`c0_summary.py` and `src_lines.py` summarize the counters and the ncu
source pages). The counters ran over the whole run. ncu profiled the
two collision launches of one step: step 8000 and step 30000 of HERMES
at 40000, step 1000 of SOLEDGE3X and of Gkeyll.

Counters, per warp of the collision kernel (32 lanes):

| Case | Species | Lanes with a slot | Charged | Busy in the sub-steps | Sub-step efficiency | Draws per sample | Warp draw passes per sub-step |
|---|---|---|---|---|---|---|---|
| HERMES 40000 | ions | 16.7 | 9.8 | 10.7 | 1.00 | 25.7 | 76.6 |
| | electrons | | | 10.6 | 0.99 | 4.9 | 13.8 |
| HERMES 160000 | ions | 26.7 | 15.8 | 16.5 | 1.00 | 25.7 | 86.3 |
| | electrons | | | 16.1 | 0.98 | 4.9 | 15.6 |
| SOLEDGE3X `self` | ions | 1.8 | 1.2 | 1.3 | 0.97 | 1 (no skew) | 1 |
| | electrons | | | 1.0 | 0.80 | 1 | 1 |
| Gkeyll | ions | 26.3 | 26.3 | 26.3 | 1.00 | 1 (no skew) | 1 |
| | electrons | | | 8.6 | 0.33 | 1 | 1 |

Sub-step efficiency is the mean `n_sub` over the mean of each warp's
largest. The other branches don't matter: no warp mixed zones with and
without a zShift table, and the A(s) branches are cheap.

ncu, share of the kernel's stall samples and busy lanes per instruction,
by source file:

| Case | `pcg32.h` (the sampler's normals) | Interpolation (`interp_hd.h`, `zones_hd.h`, `interpolate.cuh`) | Nanbu (`collisions.cu`) | Lanes, whole kernel |
|---|---|---|---|---|
| HERMES, step 8000 | 60% at 2.7 lanes | 27% at 11 to 13 | 5% at 13 | 5.2 |
| HERMES, step 30000 | 25% at 1 | 59% at 1 | 10% at 1.2 | 1.1 |
| SOLEDGE3X | 22% at 1 | 17% at 1 | 47% at 1 | 1.2 |
| Gkeyll | 25% at 7 | 34% at 32 | 35% at 12 | 11.9 |

What this shows:

- **HERMES: the skewed sampler.** Every charged particle there is
  sampled with a skew, and for the ions |eps| is large: 53% of the ion
  samples have it capped at 1 (the run's log). The bound M = 1 + |eps|
  `shape.bound` is taken from the worst case up to 5 v_th, where
  `shape.bound` is 32.4 for the ion shape and 58.8 for the electron
  shape. So an ion takes 26 draws on average, and a warp waits for
  about 80. The comment in `sample_bkg` that M is at most about 2 is
  wrong. The normals for these draws are 60% of the kernel's time, at
  under 3 busy lanes.
- **Gkeyll: the electron sub-cycles.** The electrons have 4 sub-steps
  on average, but each warp's largest is 13, so a third of the lanes are
  busy in the loop. Sorting by `n_sub` (C3) is for this.
- **The neutrals and the interpolation.** In HERMES 40% of the lanes
  with a slot are neutral, and the interpolation, which runs at the
  charged lanes' 11 to 13, is a quarter of the time. In Gkeyll it is a
  third, with every lane busy. Doing it once for both species (C2) saves
  about half of that.
- **Few particles left** (HERMES late, SOLEDGE3X). One lane per warp,
  as the spread list intends (L4b): these steps are latency, and
  nothing in this plan changes them.
- **Sub-cycles in HERMES**, the leading guess before C0: none. The ions
  always have 1 and the electrons at most 3.

So C4 goes first, and it changes how the sampler draws rather than
only its loop: with a tighter envelope, an ion sample takes about 3
draws at |eps| = 1 instead of 33 (see "A tighter envelope" below).
Then C2, then C3, which matters only for runs like the Gkeyll case.

**C1. References.**

- Seeded reference outputs from `main` for the end-to-end cases and the
  `_gpu` regression variants, as in the live-slot plan's L0, to compare
  C2 and C3 against.

Done (2026-11-27), outputs in `$SCRATCH/flan_gpu/coll/c1/refs`. Release
builds of `main` and of each phase, made after the build types were
fixed, are under `$SCRATCH/flan_gpu/coll/rel`, and each phase is
compared with a Release build of the one before it, since the build
type alone changes the last bits.

**C2. Prepare kernel.**

- A new kernel that interpolates once per step for both species and
  writes the stored values and `n_sub`. The collision kernel reads them.
  `collision_gpu`'s callers don't change.
- Tests: `test_collisions_cuda` gets a case where the split kernels and
  a run of the old path give identical slots, including with the thermal
  force on and with sub-cycling.
- Checks: `ctest` on both builds; seeded runs equal to C1's (counts and
  trajectories exactly, sums to rounding); `_gpu` regression variants
  pass their checkers; registers per thread for the new kernels from
  `-Xptxas -v` or ncu.

Done (2026-11-27), commit 022acf0, job in `$SCRATCH/flan_gpu/coll/c2`.
`collision_gpu`'s callers did change in the end: the transport loop now
calls `collisions_gpu` once per step for both species, and
`collision_gpu`, one species with its own prepare, stays for the tests.

- Registers per thread: the collision kernel went from 173 to 104, the
  prepare kernel uses about 160.
- Against C4 (Release): every count, fate and printed summary of the
  seeded runs is the same. Velocity sums differ by up to 1e-11
  relative, because the split changes where the compiler fuses
  multiply-adds and the last bits grow over the steps.
- Following impurities, HERMES `all`: 44.05 to 43.36 s at 40000
  particles, 58.41 to 55.89 s at 160000.

**C3. Collision list ordered by work** (if C0 shows the sub-cycles
matter).

- The sort after the prepare kernel, the collision kernels over the
  ordered list, neutrals skipped by the device-side count.
- Tests: a `*_cuda` test that the ordered list holds each charged live
  slot exactly once, sorted by key and in slot order within a key, with
  the neutrals counted; and that the collisions give identical slots
  with and without the order.
- Checks: as C2.

Tried and dropped (2026-11-27). The list was sorted with
`cub::DeviceRadixSort` by a key of the electron and ion sub-step buckets
(floor(log2 n_sub) + 1, 4 bits each), only when the list was packed 32
to a warp, and the collision kernels skipped the neutral and dead
entries at its start. Results were the same as C2's to rounding, and
nsys showed no gain:

| Case | Steps sorted | Collision kernel, C2 to C3 | Sort | Net |
|---|---|---|---|---|
| Gkeyll, 100000 | 1196 of 2450 | 0.654 to 0.671 s | 0.031 s | 0.05 s slower |
| HERMES `all`, 160000 | 20001 of 63121 | 7.34 to 7.21 s | 0.52 s | 0.38 s slower |

The electron sub-cycles left in a warp cost little next to the loads
every lane makes, so grouping them saves less than the sort costs. The
code is kept as a patch in `$SCRATCH/flan_gpu/coll/c3/c3.patch`.

**C4. Sampler: the tighter envelope.** Done first, after C1 (C0's
results).

- `sample_bkg` (CPU and GPU, shared in `collisions_hd.h`) draws from
  the mixture envelope and accepts with max(F, 0) / G. `skew_bound` and
  `SkewShape::bound` go, and so does the `above_bound` count, which can
  no longer happen (kept in the output as 0, or removed with the
  output's readers, to decide when writing it). The wrong comment about
  M goes with it. The `docs-site` page on the thermal force describes
  the sampling.
- Tests: in `test_thermal_force_cpu`, samples at several eps for both
  shapes reproduce the target's flow (0), heat flux and friction moment
  within their statistical error, and the mean draws per sample match
  the mixture's weight.
- Checks: `ctest` on both builds; every variant of `reg_friction` and
  `reg_thermal_force`, CPU and `_gpu`, through `check_collisions.py`;
  seeded end-to-end runs within Monte Carlo error (3 sigma) of C1's.

Done (2026-11-27), commit 66efd61, jobs in `$SCRATCH/flan_gpu/coll/rel`
(`reg*.sbatch`) and `$SCRATCH/flan_gpu/coll/c4`. `above_bound` was
removed from the output.

- Every `reg_friction` and `reg_thermal_force` variant, CPU and `_gpu`,
  passes `check_collisions.py` (39 of 41 checks, the 2 misses expected).
  The GPU z10 drift is -1206 m/s against the model's -1195; `main`
  gave -1183.
- Runs without a skew are unchanged to rounding. HERMES `all` at 4000:
  draws per ion sample 26 to 2.5, per electron sample 4.9 to 1.1;
  candidates with negative F 0.013% to 0.04% (the share quadrature of
  the clipped distribution gives, checked in `test_thermal_force_cpu`); Following impurities 43.9
  to 31.6 s. Physics summaries within Monte Carlo error.

**C5. Measure.**

- The C0 cases, `main` and the branch alternating on one node: wall time,
  nsys kernel times and ncu busy lanes. Record the results here and in
  the GPU plan's "Later". Check that runs with few particles are no
  slower (the sort's cost).

**C6. Docs and merge.**

- `code-structure.md`: the prepare kernel and the ordered list, for
  new kernels that want the stored values. Then merge to `main`.

## Decisions

All five were accepted as recommended (2026-11-27).

1. **Bit-identical results for C2 and C3.** Recommended: required, apart
   from the rounding of `atomicAdd` sums and any change in how the
   compiler fuses arithmetic (checked and stated if it happens). Only C4
   may change the random numbers.
2. **What the prepare kernel stores.** Recommended: the interpolated
   plasma values (11 doubles) and `n_sub`, rather than everything
   derived from them (`s_coeff`, eps, b for both species, about 18
   doubles). That is half the memory, and the derived values cost only
   arithmetic to recompute.
3. **When to sort.** Recommended: only when the live list is packed 32
   to a warp (`per_warp == 32`). When it is spread over the warps, warps
   have few lanes and the sort's cost isn't repaid.
4. **Sort key.** Recommended: decided by C0. If the ions' and electrons'
   `n_sub` are mostly 1 for one species, sort by the other's bucket only,
   with one list for both launches. Otherwise use one list per species.
5. **The CPU.** Recommended: unchanged apart from C4's shared sampler.
   A CPU thread doesn't wait for other lanes, so the split and the sort
   don't help it.
