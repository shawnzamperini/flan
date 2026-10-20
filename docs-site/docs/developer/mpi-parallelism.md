# MPI Parallelism

Flan runs on several MPI processes (*ranks*) by giving each one a share of the
impurity particles to follow. The particles don't interact with each other,
so the ranks work independently for almost the whole run. They only
communicate twice: once at the start, to share the background plasma, and
once at the end, to add up their results. Within each rank, the particles
are followed in parallel with OpenMP threads on the CPU (see
[OpenMP Layout](openmp-layout.md)) or with the GPUs the rank can see.

## Overview

```
rank 0                         rank 1 ... rank N-1
------                         -------------------
load options                   load options
read background
broadcast background   ----->  receive background
follow imp_num/N particles     follow imp_num/N particles
reduce statistics      <-----  send statistics
calculate derived quantities
save <case_name>.nc
```

Every rank calls `flan()`, so every step below happens on every rank unless
it says otherwise.

## Startup

**Options.** Each rank loads the options itself (`Options::load_options`).
The input options are compiled into the input file, so every rank already
has them and nothing needs to be sent.

**Background.** Only rank 0 reads the background (`Background::read_bkg`).
Reading a Gkeyll background calls postgkyl through `read_gkyl.py` and writes
intermediate files, which only needs to happen once. Rank 0 then sends it to
the other ranks with `Background::broadcast`. That broadcasts the dimensions,
the times, cell centers and cell edges, and every 3D and 4D array, using
`MPI_Bcast` from rank 0 (`Utilities::mpi_broadcast_vector` and
`Vector3D/Vector4D::broadcast`). Every rank ends up with a full copy of the
background, so the memory used per rank doesn't shrink with more ranks.

**OpenADAS.** If ionization/recombination is on, every rank reads the
OpenADAS rate files itself, one rank at a time, with an `MPI_Barrier` between
them so they don't all read the same files at once.

## Following the particles

`ImpurityTransport::main_loop` divides `imp_num` between the ranks:

```cpp
int parts_per_rank {opts.imp_num() / nprocs};

// Give last process the remainder.
if (rank == nprocs - 1) parts_per_rank += opts.imp_num() % nprocs;
```

With GPUs on, each rank divides its share again between the GPUs it can see,
with GPU 0 getting the remainder. Each rank then runs the main loop on its
particles with no communication. Ranks can finish at different times. They
wait for each other at the reduction afterwards.

Each rank needs its own random numbers, or every rank would follow the same
particles. The random number generator seeds include the rank:

| | Seed | Stream |
|---|---|---|
| CPU | `seed + rank` | OpenMP thread number |
| GPU | `seed + rank * num_gpus + gpu` | slot index |

With the default `seed = -1` (any negative value), the clock is used instead
of `seed`, with the rank still added, so every run is different.

## Combining the results

After the main loop, every rank has its own `Impurity::Statistics`: counts,
weights and weighted sums on the background grid. These are summed onto
rank 0 with `Impurity::reduce_stats`:

1. Each rank packs its ten statistics arrays (counts, weights, charge, `s`
   and the six velocity components) into one buffer of doubles.
2. `MPI_Reduce` with `MPI_SUM` adds the buffers onto rank 0.
3. Rank 0 unpacks the sums into a new `Statistics` object. The other ranks
   get an empty one.

The statistics are sums over particles, so summing them across ranks gives
exactly the same result as following every particle on one rank. Only
rank 0 then calculates the density, average velocities, average charge and
average `s` from them, saves the output file and prints the timing summary.

## Running with MPI

Launch Flan with your system's MPI launcher, e.g.:

```bash
OMP_NUM_THREADS=5 mpiexec -np 2 ./test_exb
```

or `srun` on a cluster (see [Running Flan](../running/running-flan.md) and
[Perlmutter](../running/running-perlmutter.md)). Each rank uses
`OMP_NUM_THREADS` OpenMP threads, so the total number of cores used is the
number of ranks times `OMP_NUM_THREADS`.

!!! note "GPUs and ranks"
    Each rank uses every GPU it can see (`cudaGetDeviceCount`). With more than
    one rank on a node, limit each rank to its own GPUs (e.g., with
    `CUDA_VISIBLE_DEVICES` or your scheduler's GPU binding), or the ranks will
    share the same GPUs. One rank per node, using all of that node's GPUs, is
    the simplest setup.

## Testing

`tests/test_mpi.cpp` checks the MPI communication: broadcasting vectors, the
`Vector3D`/`Vector4D` arrays and the full background (including the
simulation bounds), and summing the statistics with `reduce_stats`. `ctest`
runs it on 3 ranks. To run it by hand:

```bash
mpiexec -n 3 ./tests/test_mpi
```

Only rank 0 prints the full GoogleTest output. The other ranks print just
their failures, tagged with the rank, e.g. `[rank 1] ...: Failure`. A change
to what's broadcast or reduced should come with a check here.

## Things to know

- **Output.** Only rank 0 prints the startup banner, the progress updates
  and the timing summary. The progress updates only count rank 0's
  particles, which is a good estimate of the whole run since every rank has
  about the same number.
- **Warnings.** The end-of-run warnings (ionization/recombination
  probability above 1, particles past `imp_max_time`) are printed by each
  rank for its own particles, not summed across ranks.
- **Errors.** An error anywhere stops every rank. `Errors::fatal_error`
  calls `MPI_Abort`, so the other ranks don't hang waiting at the next
  broadcast or reduction.
- **MPI setup.** The input file must call `MPI_Init` before `flan()` and
  `MPI_Finalize` after it.
