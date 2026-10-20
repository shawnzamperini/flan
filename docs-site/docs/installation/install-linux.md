# Linux Installation

These are the steps for installing Flan on a Linux machine, such as a laptop
or workstation. See the [Installation](installation.md) overview for more on
each step and on choosing between a CPU and GPU build.

## 1. Install the requirements

You need a C++23 compiler, CMake 3.26 or newer, MPI, git and conda. On
Ubuntu, for example:

```bash
sudo apt install build-essential cmake git openmpi-bin libopenmpi-dev
```

If your distribution's CMake is older than 3.26, a newer one can be installed
with `pip install cmake` or from [cmake.org](https://cmake.org/download/).
For conda, [Miniconda](https://docs.conda.io/en/latest/miniconda.html) is
enough.

For a GPU build, you also need an NVIDIA GPU with its driver, and the
[CUDA Toolkit](https://developer.nvidia.com/cuda-downloads). Check that both
work with:

```bash
nvidia-smi
nvcc --version
```

## 2. Clone Flan

```bash
git clone https://github.com/shawnzamperini/flan.git
cd flan
```

## 3. Build the dependencies

From the top-level `flan/` directory, run:

```bash
machines/mkdeps.linux.sh
```

This downloads and builds zlib, HDF5 and NetCDF in `$HOME/flansoft`, creates
the `flan` conda environment with postgkyl, and writes the file that tells
CMake where everything is. It takes a while. To put the dependencies
somewhere else, change this line in `machines/mkdeps.linux.sh` before running
it:

```bash
export FLANSOFT=$HOME/flansoft
```

Then activate the environment. Do this before building, and every time
before running Flan:

```bash
conda activate flan
```

## 4. Build, install and test Flan

From the top-level `flan/` directory, with the `flan` environment active:

For a CPU build:

```bash
mkdir build_cpu && cd build_cpu
cmake ..
make -j8
make install
ctest
```

For a GPU build:

```bash
mkdir build_gpu && cd build_gpu
cmake .. -DENABLE_CUDA=ON
make -j8
make install
ctest
```

`make` also creates `flan/scripts/flan_setup.sh`, and `make install` puts
`libflan.so` in `flan/lib/`. `ctest` should end with
`100% tests passed`. If your GPU isn't an A100 or another Ampere GPU (compute
capability 8.0 or 8.6), see
[Supported GPUs](installation.md#cpu-or-gpu-build) before building.

You're ready to set up a case, see [Running Flan](../running/running-flan.md).

## Troubleshooting

See [Troubleshooting](installation.md#troubleshooting) for common problems.
For anything else, contact Shawn Zamperini
([zamperinis@fusion.gat.com](mailto:zamperinis@fusion.gat.com)).
