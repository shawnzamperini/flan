# Perlmutter Installation

These are the steps for installing Flan on NERSC's Perlmutter. See the
[Installation](installation.md) overview for more on each step.

Perlmutter's GPU nodes have NVIDIA A100 GPUs, which the GPU build supports,
so the GPU build is the one to use. It runs on Perlmutter's CPU nodes too,
with `use_gpu = "off"`. On Perlmutter, Flan is built with the Cray compiler
wrappers (`cc` and `CC`), which come with Cray MPICH.

## 1. Clone Flan

```bash
git clone https://github.com/shawnzamperini/flan.git
cd flan
```

## 2. Build the dependencies

Load the conda module, then run the Perlmutter dependency script from the
top-level `flan/` directory:

```bash
module load conda
machines/mkdeps.perlmutter.sh
```

This builds zlib, HDF5 and NetCDF with the Cray wrappers in `$HOME/flansoft`,
creates the `flan` conda environment with postgkyl, and writes the file that
tells CMake where everything is. To put the dependencies somewhere else,
change this line in `machines/mkdeps.perlmutter.sh` before running it:

```bash
export FLANSOFT=$HOME/flansoft
```

Then activate the environment. Do this before building, and every time
before running Flan, including in job scripts:

```bash
module load conda
conda activate flan
```

## 3. Build, install and test Flan

From the top-level `flan/` directory, with the `flan` environment active:

```bash
mkdir build_gpu && cd build_gpu
cmake .. -DENABLE_CUDA=ON
make -j8
make install
```

Check that the CMake output says `NVIDIA GPU support enabled`. `make` also
creates `flan/scripts/flan_setup.sh`, and `make install` puts `libflan.so` in
`flan/lib/`.

Run the tests on a GPU node, e.g., in an interactive job, since they need a
GPU and run on several MPI ranks:

```bash
ctest
```

It should end with `100% tests passed`.

For a CPU-only build, e.g., if you'll only use the CPU nodes, use
`mkdir build_cpu && cd build_cpu` and `cmake ..` instead, then the same
`make`, `make install` and `ctest`.

You're ready to set up a case, see [Running Flan](../running/running-flan.md)
and [Running on Perlmutter](../running/running-perlmutter.md).

## Troubleshooting

See [Troubleshooting](installation.md#troubleshooting) for common problems.
For anything else, contact Shawn Zamperini
([zamperinis@fusion.gat.com](mailto:zamperinis@fusion.gat.com)).
