# Running Flan

## Summary

Setting up a Flan case from scratch consists of the following steps:

- Create case directory: `flan_setup.sh`
- Modify input file
- Compile input file with Flan
- Run case/submit job 

## Setup

A built-in test case that demonstrates the ExB drift with a single particle
will be used to show how to run Flan. First create the run directory with the 
command:

```bash
/path/to/flan/scripts/flan_setup.sh test_exb
```

This will create a directory called test\_exb and populate it with an empty
input file as well as the needed files to compile it with Flan. We will not
cover the various input options here, and will instead copy over the input file
used to run the ExB drift test:

```bash
cd test_exb
cp /path/to/flan/tests/test_exb.cpp .
```

The input file must have the same name as the directory because it is baked 
into the Makefile that was generated to compile the input file. Flan must be
run in the conda environment that was created as part of the installation 
process because it contains some needed paths, primarily for interfacing with 
Gkeyll data.  Ensure it's activated with `conda activate flan`. You typically 
will see `(flan)` at the start of your command line when it's activated. Next 
make the input file with just `make`.

Running the simulation depends on if you can run it from the command line or
need to use a job scheduler. Running from the command line will be covered 
first, and job submission scripts will be provided after. 

# Running from Command Line

To run from the command line simply run the executable:

```bash
./test_exb
```

By default, this will run with one MPI process and likely the maximum number of
OpenMP threads available. On smaller machines like laptops or workstations only
one MPI processes is needed since it is mainly for distributed computing across
multiple compute nodes, such as on clusters. To explicitly set the number of
OpenMP threads used by Flan, e.g., to 5, use:

```bash
OMP_NUM_THREADS=5 ./test_exb
```

You may often need to reduce the number of OpenMP threads to stay within memory
limits of your machine, as each thread comes with memory overhead.

To control the number of MPI processes, e.g., to use 2, use:

```bash
mpiexec -np 2 ./test_exb
```

Alternative commands include `mpirun` and `srun`, depending on your system.

Control the number of OpenMP threads available to each process the same as
before:

```bash
OMP_NUM_THREADS=5 mpiexec -np 2 ./test_exb
```

# Running with Job Submission Script

For convenience job submission scripts are provided for the following machines:

- [Perlmutter](running-perlmutter.md)
- [Omega](running-omega.md)

