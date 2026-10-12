# Perlmutter Installation

Installing Flan on Perlmutter involves downloading and compiling
the dependencies and then compiling the code via CMake. It uses the NERSC 
provided Cray C/C++ compiler wrappers (`cc`/`CC`).

## Install Dependencies

Navigate to the top-level `flan/` directory and run:

```bash
machines/mkdeps.perlmutter.sh
```

This will download and install all the dependencies in `$HOME/flansoft` by 
default. If you wish to have them installed elsewhere, change the following
line in `mkdeps.perlmutter.sh` to the desired directory:

```bash
export FLANSOFT=/path/to/dir/flansoft
```

This script also handles creating the conda environment that Flan must run in.
You must activate the environment before running Flan by first loading the
conda module with `module load conda` and then activating the environment with
`conda activate flan`. 

## Install Flan

Navigate to `flan/build`. This is where Flan will be compiled as a library. 
Run CMake to create the Makefile, then make and install the library.

```bash
cd build
cmake ..
make && make install
```

Flan will be installed at `flan/lib/libflan.so`. 

### Troubleshooting

Please reach out to a developer if you encounter issues so the instructions 
can be updated.
