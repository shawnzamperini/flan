# Running on Perlmutter

On Perlmutter, compile the case on a login node (see
[Running Flan](running-flan.md)), then submit the run as a Slurm job. Below
are example job scripts for the CPU nodes and the GPU nodes. Change
`--account` to your NERSC project, and the executable name to your case.
Copies of both are in the repository, as
`scripts/submit-flan-perlmutter-cpu.sbatch` and
`scripts/submit-flan-perlmutter-gpu.sbatch`.

!!! note "Starting points"
    These scripts are a place to start, not tuned settings. The CPU script is
    the one Flan has been run with. The GPU script hasn't been tested yet.
    How many CPUs to ask for and how to bind the threads depends on the job,
    see NERSC's guide to
    [process and thread affinity](https://docs.nersc.gov/jobs/affinity/).

## CPU nodes

This script runs one MPI rank per node with 256 OpenMP threads each, on 2
nodes. `--cpus-per-task` counts hardware threads, and a CPU node has 128 cores
with 2 hardware threads each, so 256 gives the rank the whole node:

```bash
#!/bin/bash -l
#SBATCH --job-name=my_case
#SBATCH --account=mXXXX         # Your NERSC project
#SBATCH --constraint=cpu
#SBATCH --nodes=2
#SBATCH --ntasks-per-node=1     # One MPI rank per node
#SBATCH --cpus-per-task=256     # The whole node (128 cores x 2 threads)
#SBATCH --qos=regular           # debug (one node, short runs) or regular
#SBATCH --time=01:00:00

# OpenMP settings
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
export OMP_PLACES=cores         # Bind threads to cores
export OMP_PROC_BIND=spread     # Spread threads across cores

# Threads available to Python (postgkyl, when reading a Gkeyll background)
export NUMEXPR_MAX_THREADS=$SLURM_CPUS_PER_TASK

# The flan environment is needed on the compute node too
module load conda
conda activate flan

srun --cpu-bind=cores ./my_case
```

## GPU nodes

Each GPU node has 4 NVIDIA A100 GPUs. Flan uses every GPU a rank can see, so
run one rank per node and give it all 4. Set `inpts["use_gpu"] = "cuda"` in
the input file, and compile against a GPU build of Flan
(see [Perlmutter Installation](../installation/install-perlmutter.md)).

```bash
#!/bin/bash -l
#SBATCH --job-name=my_case
#SBATCH --account=mXXXX         # Your NERSC project
#SBATCH --constraint=gpu
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1     # One MPI rank per node, using all 4 GPUs
#SBATCH --gpus-per-node=4
#SBATCH --cpus-per-task=128     # The whole node (64 cores x 2 threads)
#SBATCH --qos=regular
#SBATCH --time=01:00:00

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
export OMP_PLACES=cores
export OMP_PROC_BIND=spread
export NUMEXPR_MAX_THREADS=$SLURM_CPUS_PER_TASK

module load conda
conda activate flan

srun --cpu-bind=cores ./my_case
```

A GPU node has 64 cores, so 128 hardware threads is the whole node. The
OpenMP threads are still used on the GPU nodes, for reading the background
and other setup.

## Submitting

If the script is saved as `submit.sh`:

```bash
sbatch submit.sh        # Submit the job
squeue -u $USER         # Check on it
```

Flan's output goes to `slurm-<jobid>.out` in the directory the job was
submitted from, and the results to `<case_name>.nc` in the case directory.
Use the `debug` QOS for quick tests: it's limited to short, small jobs, but
usually starts much sooner.

See [Running Flan](running-flan.md#troubleshooting) for troubleshooting.
