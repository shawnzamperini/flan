# Running on Omega

Omega is General Atomics' cluster. Flan runs on Omega's CPU nodes, set up and
compiled the same as on any Linux machine, with a few differences.

## Installing on Omega

Follow the [Linux](../installation/install-linux.md) instructions, but build
the dependencies with the Omega script:

```bash
machines/mkdeps.omega.sh
```

It puts the dependencies in `/fusion/projects/codes/flan/$USER/flansoft`
instead of your home directory, and it creates the `flan` environment with
**mamba** instead of conda, since Omega prefers mamba. Activate it with:

```bash
mamba activate flan
```

Omega has several compiler modules, and others loaded can get in the way.
Load GCC 13 before building, the same as in the job script below:

```bash
module purge
module load gcc/13.x
```

## Job script

Compile the case on a login node (see [Running Flan](running-flan.md)), then
submit the run as a Slurm job. This script runs one MPI rank with 16 OpenMP
threads:

```bash
#!/bin/bash -l
#SBATCH --job-name=my_case
#SBATCH --export=ALL
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem-per-cpu=15G
#SBATCH --partition=short
#SBATCH --time=00:30:00

# Other modules get in the way, purge them
module purge
module load gcc/13.x

# The flan environment is needed on the compute node too
eval "$(mamba shell hook --shell bash)"
mamba activate flan

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

./my_case

# Copy the results and the saved background files to the results area
export RESULTSDIR="/fusion/projects/results/flan/$USER/my_case"
mkdir -p $RESULTSDIR
cp bkg* $RESULTSDIR
cp *.nc $RESULTSDIR
```

Change `my_case` to your case name. The last few lines copy the results
(`<case_name>.nc`) and the saved Gkeyll background files to Omega's results
area, which is where results are meant to be kept.

## Submitting

If the script is saved as `submit.sh`:

```bash
sbatch submit.sh        # Submit the job
squeue -u $USER         # Check on it
```

Flan's output goes to `slurm-<jobid>.out` in the directory the job was
submitted from.

See [Running Flan](running-flan.md#troubleshooting) for troubleshooting.
