# OpenMP Layout

On the CPU, Flan follows the particles on each MPI rank in parallel with
OpenMP threads. Each step of the main loop is a loop over the particle
slots, and OpenMP splits that loop between the threads. This page describes
how that's laid out, and how the few places where threads share something
are handled. For running on several ranks, see
[MPI Parallelism](mpi-parallelism.md). With GPUs on, the main loop runs on
the GPUs instead and OpenMP only helps with setup.

## Threads per rank

The number of threads is set with `OMP_NUM_THREADS`, the same for every
rank:

```bash
OMP_NUM_THREADS=8 ./my_case                # 1 rank, 8 threads
OMP_NUM_THREADS=8 mpiexec -n 4 ./my_case   # 4 ranks, 8 threads each
```

Flan prints the number it's using at startup (`OMP Threads/Rank`). On a
cluster, bind the threads to cores so they don't move around, e.g., with
`OMP_PLACES=cores` and `OMP_PROC_BIND=spread` (see
[Perlmutter](../running/running-perlmutter.md)).

## The main loop

Every CPU version of a main loop step (`*_cpu`, see
[Code Structure](code-structure.md#the-main-loop)) has the same shape: a
loop over the slots, split between the threads, that skips dead slots and
updates each alive particle on its own:

```cpp
#pragma omp parallel for
for (int i = 0; i < slots.N(); ++i)
{
	// Skip dead particles
	if (slots.state()[i] > 0) continue;

	// Update particle i
	...
}
```

The particles don't interact with each other, so each thread only ever reads
and writes its own slots, and no locking is needed. The default `static`
schedule splits the slots into one contiguous block per thread, which keeps
each thread's memory accesses together. The slots are stored as a structure
of arrays (one array per quantity), so a thread walks through each array in
order.

The parallel region starts and ends inside each step, so the threads wait
for each other at the end of every step before the next one starts. The
steps run in a fixed order and each depends on the one before (e.g., the
Boris push needs the cells found in the previous step), so this is needed.

## Random numbers

Collisions, ionization/recombination, the core boundary condition and
launching new particles all need random numbers. Each thread has its own
PCG32 generator, so threads never share one:

```cpp
#pragma omp parallel
{
	// Grab our RNG for this thread, each thread has its own (and they're
	// seeded uniquely).
	int tid = omp_get_thread_num();
	pcg32& rng = rngs[tid];

	#pragma omp for
	for (int i = 0; i < slots.N(); ++i)
	{
		...
	}
}
```

`follow_impurities` creates one generator per thread before the main loop
(`omp_get_max_threads()` of them). They all use the same seed (`seed` plus
the MPI rank) but a different stream (the thread number), which gives every
thread its own independent sequence. Using `#pragma omp parallel` and then
`#pragma omp for`, instead of `#pragma omp parallel for`, is what lets each
thread look up its generator once, before the loop starts.

## Shared results

A few things are added up across all the particles. These are the only
places threads write to the same memory.

| What | Where | How |
|---|---|---|
| Statistics (counts, weights, velocities, charge, `s`) | `Impurity::record_stats_cpu` | `#pragma omp atomic update` on each add |
| Ionization/recombination probability > 1 warnings | `OpenADAS::ioniz_recomb` | `reduction(+:...)` |
| Particles past `imp_max_time` | `step_cpu` | `#pragma omp atomic` |
| Dead slots | `Slots::fill_slots_cpu` | `reduction(+:dead_count)` |
| Slots claimed for new particles | `Slots::fill_slots_cpu` | `#pragma omp atomic capture` |
| Error messages (e.g., `s = nan`) | `Collisions::nanbu_coll` | `#pragma omp critical` around the printing |

**Statistics.** Many particles can be in the same cell, so each add to the
statistics arrays is atomic. Atomics are cheap when threads rarely add to the
same cell at the same time, which is the usual case since there are many more
cells than threads. The commented-out *buffered* version in
`impurity_stats.cpp` gives each thread its own copy of the arrays and adds
them up at the end instead. It was slower in testing, but could win when
there are far more particles than cells.

**Filling slots.** `fill_slots_cpu` first counts the dead slots, and decides
how many to revive (all of them, or the number of particles left to launch,
whichever is smaller). Then each thread claims dead slots by incrementing a
shared counter with `atomic capture`, and stops reviving once the counter
reaches that number. That way exactly the right number of particles is
launched, without a lock.

## Setup

OpenMP is also used outside the main loop, mostly in reading the background:
converting and checking the Gkeyll data in `read_gkyl.cpp`, filling the test
backgrounds in `read_test.cpp`, and calculating the derived quantities in
`impurity_stats.cpp`. These are simple `parallel for` loops over grid cells.

## Things to know

- **No MPI inside parallel regions.** Flan calls `MPI_Init`, not
  `MPI_Init_thread`, so MPI may only be called from outside OpenMP parallel
  regions. This is also why the `nan` checks in the collisions only print a
  warning instead of calling `Errors::fatal_error`, which calls `MPI_Abort`.
- **Reproducibility.** With a fixed `seed`, a run isn't exactly repeatable
  with threads. Which thread (and so which random number sequence) handles
  which particle depends on the number of threads, which slots get new
  particles depends on the order threads reach the shared counter, and the
  atomic adds to the statistics happen in a different order each run, which
  changes the last few digits. The results are statistically the same. For a
  repeatable run, use one thread and one rank.
- **Adding a step.** A new CPU step should follow the same pattern: a
  `parallel for` over the slots that skips dead ones, the thread's own
  generator if it needs random numbers, and an atomic, reduction or critical
  section for anything shared between particles.
