# Linux Installation

Installing Flan on a general Linux machine involves downloading and compiling
the dependencies and then compiling the code via CMake. 

## Install Dependencies

Navigate to the top-level `flan/` directory and run:

```bash
machines/mkdeps.linux.sh
```

This will download and install all the dependencies in `$HOME/flansoft` by 
default. If you wish to have them installed elsewhere, change the following
line in `mkdeps.linux.sh` to the desired directory:

```bash
export FLANSOFT=/path/to/dir/flansoft
```

This script also handles creating the conda environment that Flan must run in.
You must activate the environment before running Flan with 
`conda activate flan`.

## Install Flan

Navigate to `flan/build`. This is where Flan will be compiled as a library. 
Run CMake to create the Makefile, then make and install the library.

```bash
cd build
cmake ..
make && make install
```

Flan will be installed at `flan/lib/libflan.so`. This also creates a bash
script `scripts/flan_setup.sh`, which is used when setting up a new case.

### Troubleshooting

There are probably numerous issues that can pop up, as each machine is 
different. Please reach out to a developer if you encounter issues so the 
instructions can be updated.

