#!/bin/bash -l
#SBATCH --job-name=flan-case
#SBATCH --account=m3739
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --constraint=cpu
#SBATCH --cpus-per-task=128
#SBATCH --partition=debug
#SBATCH --time=00:30:00

# OpenMP settings
export OMP_NUM_THREADS=128      # Match --cpus-per-task above
export OMP_PLACES=cores			# Bind threads to cores
export OMP_PROC_BIND=spread     # Spread threads across cores

# This controls the number of threads available in python. Make it match
# cpus-per-task
export NUMEXPR_MAX_THREADS=128

# Initialize conda for compute node
module load conda
conda activate flan

# Run Flan run
srun --cpu-bind=cores ./flan-case
