# Perlmutter Batch Script

The following is an example batch script used to submit a Flan job on 
Perlmutter. 

```bash
#!/bin/bash -l
#SBATCH --job-name=my-job-name
#SBATCH --account=m3739     # Update with project number
#SBATCH --nodes=2           # 2 nodes
#SBATCH --ntasks=2          # 2 MPI tasks
#SBATCH --constraint=cpu
#SBATCH --cpus-per-task=128 # Use all CPUs on node
#SBATCH --qos=regular       # debug (only one node allowed) or regular
#SBATCH --time=01:00:00     # Run for an hour

# OpenMP settings
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK      # Match --cpus-per-task above
export OMP_PLACES=cores			# Bind threads to cores
export OMP_PROC_BIND=spread     # Spread threads across cores

# This controls the number of threads available in python. Make it match
# cpus-per-task
export NUMEXPR_MAX_THREADS=$SLURM_CPUS_PER_TASK

# Initialize conda for compute node
module load conda
conda activate flan

# Run Flan run
srun --cpu-bind=cores ./my_flan_case
```
Compile the input file as normal, using the Makefile provided when setting up 
a case directory. If the batch script is saved as `submit-flan-perlmutter.sh`, 
then submit the job with:

```bash
sbatch submit-flan-perlmutter.sh
```

Output is put into a generically named slurm file. 
