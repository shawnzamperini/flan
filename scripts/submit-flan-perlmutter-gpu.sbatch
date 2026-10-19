#!/bin/bash -l
#SBATCH --job-name=flan-case
#SBATCH --account=m3739
#SBATCH --nodes=1
#SBATCH --ntasks=1              # One MPI rank per node, using all 4 GPUs
#SBATCH --constraint=gpu
#SBATCH --gpus-per-node=4
#SBATCH --cpus-per-task=128     # The whole node (64 cores x 2 threads)
#SBATCH --partition=debug
#SBATCH --time=00:30:00

# Flan must be built with -DENABLE_CUDA=ON, and the input file must set
# inpts["use_gpu"] = "cuda", or the run stays on the CPU.

# OpenMP settings. The threads are still used on the GPU nodes for reading
# the background and other setup.
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
