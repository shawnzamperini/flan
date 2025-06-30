#!/bin/bash -l
#SBATCH --job-name=flan-test
#SBATCH --export=ALL
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem-per-cpu=15G
#SBATCH --nodes=1
#SBATCH --time=00:30:00
#SBATCH --partition=short

# Other modules get in the way, purge them
module purge
module load gcc/13.x

# Run Flan run
./test

# Be a good user: Move results to the results directory
export RESULTSDIR="/fusion/projects/results/flan/zamperinis/test"
mkdir $RESULTSDIR
cp bkg* $RESULTSDIR
cp *.nc $RESULTSDIR

