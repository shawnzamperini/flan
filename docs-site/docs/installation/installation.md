# Installation

Flan is a C++ library, `libflan.so`, built with CMake. Installing it takes
three steps:

1. Clone the repository.
2. Build the dependencies and create the `flan` conda environment, with the
   script for your machine.
3. Configure, build, install and test Flan with CMake.

The steps are the same everywhere, apart from the dependency script and a few
machine-specific settings. The [Linux](install-linux.md) and
[Perlmutter](install-perlmutter.md) pages have the exact commands for each.
This page explains the choices and what each step does.

## CPU or GPU build

Flan can be built two ways:

- **CPU build.** Particles are followed in parallel with MPI and OpenMP.
  This works on any Linux machine.
- **GPU build.** Everything in the CPU build, plus CUDA versions of the main
  loop that follow the particles on GPUs. This needs an **NVIDIA GPU** and the
  CUDA Toolkit, so it doesn't work with AMD or Intel GPUs.

A GPU build can still run on the CPU: whether a simulation uses the GPUs is
chosen at run time with the `use_gpu` input option (see
[Input Options](../input-options.md)). So if your machine has an NVIDIA GPU,
the GPU build is the one to use.

!!! note "Supported GPUs"
    The GPU build is compiled for compute capability 8.0 and 8.6, which
    covers the NVIDIA A100 (Perlmutter) and other Ampere GPUs, such as the
    RTX 30 series and RTX A series. For a different GPU, change
    `CMAKE_CUDA_ARCHITECTURES` in the top-level `CMakeLists.txt` to your GPU's
    compute capability, e.g., `90` for an H100. You can find it with
    `nvidia-smi --query-gpu=name,compute_cap --format=csv`.

## Requirements

Your machine needs to provide:

- A C++ compiler that supports C++23 (GCC 11.4 has been tested)
- CMake 3.26 or newer
- MPI (e.g., OpenMPI, MPICH or Cray MPICH)
- OpenMP (included with GCC)
- git
- conda (e.g., Miniconda) or mamba
- An internet connection when building, to download the dependencies and
  GoogleTest
- **GPU build only:** an NVIDIA GPU and the CUDA Toolkit (CUDA 12.8 has been
  tested)

zlib, HDF5 and NetCDF (C and C++) are downloaded and built by the dependency
script, along with postgkyl, which Flan uses to read Gkeyll data.

## Clone the repository

```bash
git clone https://github.com/shawnzamperini/flan.git
cd flan
```

## Build the dependencies

Run the dependency script for your machine from the top-level `flan/`
directory, e.g., for a general Linux machine:

```bash
machines/mkdeps.linux.sh
```

This builds zlib, HDF5 and NetCDF in `$HOME/flansoft` (change `FLANSOFT` in
the script to put them elsewhere). It also creates the `flan` conda
environment, which has postgkyl and the Python packages Flan needs, and sets
the environment variables Flan uses to find its Python scripts. Finally, it
writes `install-deps/build-opts.cmake`, which tells CMake where the
dependencies are. CMake won't configure without it, so the dependencies
always come first.

Activate the environment before building and before running Flan:

```bash
conda activate flan
```

!!! note "Build and run in the same environment"
    The `flan` environment includes its own MPI. Activate it before running
    CMake, so Flan is built with the same MPI it'll run with. Building with
    one MPI and running with another usually crashes at the first MPI call.

## Build Flan

The CPU and GPU builds each go in their own build directory, so you can have
both. From the top-level `flan/` directory:

**1. Make a build directory and configure.** CMake checks the compilers and
dependencies and writes the Makefiles:

For a CPU build:

```bash
mkdir build_cpu && cd build_cpu
cmake ..
```

For a GPU build:

```bash
mkdir build_gpu && cd build_gpu
cmake .. -DENABLE_CUDA=ON
```

The output should end with `-- Build files have been written to: ...`. For a
GPU build, check that it says `NVIDIA GPU support enabled`.

**2. Build.** This compiles the library and the tests, and creates
`flan/scripts/flan_setup.sh`, which sets up new cases (see
[Running Flan](../running/running-flan.md)). `-j8` builds with 8 processes.
Use however many cores you have.

```bash
make -j8
```

**3. Install.**

```bash
make install
```

This copies `libflan.so` to `flan/lib/`, and the CMake files that cases use
to find Flan to `flan/lib/cmake/flan/`.

!!! note "One installed build at a time"
    Both builds install to the same `flan/lib/`, so the last `make install`
    is the one cases use. If you switch between a CPU and GPU build, run
    `make install` in the one you want, then rebuild your case.

**4. Run the tests.**

```bash
ctest
```

All the tests should pass:

```
100% tests passed, 0 tests failed out of 7
```

A CPU build has 7 tests. A GPU build has those plus the GPU tests, 17 in
all, and they need a GPU to run on. One of the tests checks the MPI
communication, so it runs Flan's test program on 3 MPI ranks. If a test
fails, `ctest --output-on-failure` shows what went wrong.

### Build options

These are set with `-D` when configuring, e.g., `cmake .. -DENABLE_CUDA=ON`.

| Option | Default | Description |
|---|---|---|
| `ENABLE_CUDA` | `OFF` | Build the GPU version. Needs the CUDA Toolkit. |
| `CMAKE_BUILD_TYPE` | `Release` | `Release` for normal use. `Debug` turns off optimization and turns on extra checks, e.g., the GPU waits after every kernel and prints any error. `RelWithDebInfo` is optimized but keeps debugging information. |

### Rebuilding

After pulling new changes, rebuild and reinstall from the build directory:

```bash
make -j8 && make install
```

If files were added to or removed from Flan, run `cmake ..` first so CMake
picks them up.

## Troubleshooting

Every machine is a little different, so problems can come up that these
instructions don't cover. A few common ones:

- **CMake can't find `build-opts.cmake`:** the dependency script hasn't been
  run yet, or it failed partway. Run it again and check its output.
- **CMake can't find the CUDA compiler:** the CUDA Toolkit isn't installed or
  `nvcc` isn't on your `PATH`.
- **A crash right at the start of a run, mentioning MPI:** Flan was probably
  built and run with different MPIs. See the note on building and running in
  the same environment above.

For anything else, contact Shawn Zamperini
([zamperinis@fusion.gat.com](mailto:zamperinis@fusion.gat.com)), so the
instructions can be updated too.

## Machine instructions

- [Linux](install-linux.md): any Linux machine, such as a laptop or
  workstation.
- [Perlmutter](install-perlmutter.md): NERSC's Perlmutter.
