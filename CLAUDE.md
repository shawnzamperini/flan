# Flan

Monte Carlo impurity transport code: follows impurity ions (mostly
tungsten) through a plasma background from Gkeyll, SOLEDGE3X, or the
built-in slab, with Boris pushes, Nanbu Coulomb collisions (including the
thermal force via heat-flux-skewed sampling), OpenADAS ionization and
recombination, and a sputtered source from RustBCA yields. C++23 with
OpenMP and MPI; optional CUDA (no thermal force, sputtering or sheath on the GPU yet).

Licensed GA noncommercial. Other codes (e.g. OpenEdge, GPLv3) may be used
for ideas only: never copy their code.

## Layout

| Path | Contents |
|---|---|
| `src/`, `include/` | The library. `src/cuda/` holds the GPU versions. |
| `tests/` | Unit tests (GoogleTest), built with the library, run with `ctest`. |
| `regression/` | Physics regression cases `reg_*.cpp` on the slab background, and the checkers `check_collisions.py`, `check_drifts.py`. See `regression/README.md`. |
| `python/` | Plotting and analysis scripts. |
| `docs-site/` | MkDocs documentation. All docs go here. |
| `install-deps/` | Dependency build scripts and `build-opts.cmake`. |

## Building and testing

- Two build directories in the repo root: `build_cpu` (CPU only) and
  `build_gpu` (`-DENABLE_CUDA=ON`, Release). Build with `make -j10`, test
  with `ctest -j8` in each. Each test runs with `FLAN_TEST_THREADS` (CMake
  cache, default 4) OpenMP threads, so `-j` doesn't oversubscribe the
  node. `MPITest` needs `srun`, so it only runs on a compute node.
- Sources are globbed (`file(GLOB ...)`): after adding a `.cpp` or `.cu`
  file, re-run `cmake .` in the build directory.
- `make install` updates `lib/libflan.so`, which installed inputs link
  against. Plain `make` does not.
- Tests named `*_cpu` always build; `*_cuda` only with CUDA. A CUDA test
  that defines its own kernels is a `.cu` file, otherwise `.cpp`.
- After changes to collisions, the push or the slab background, run the
  regression cases and their checkers (`regression/README.md`). They take
  from about a second (drifts) to several minutes (collisions).

## Conventions

- Comment code well, Doxygen style (`/** ... */` with `@param`,
  `@return`), saying why as well as what.
- An input option touches many files that don't check each other; follow
  "Adding an input option" in `docs-site/docs/developer/code-structure.md`
  (options.h/.cpp, load_options.cpp, save_results.cpp, and for GPU use
  options_device.h and cuda/options.cu), and document it in
  `docs-site/docs/input-options.md`.
- New physics gets a regression case on the slab background (set entirely
  by `slab_*` input options, nothing hard-coded) with an analytic or fluid
  comparison in a checker script.
- CPU and GPU versions of the same routine should stay in step; GPU
  changes should be checked with a `*_gpu` regression variant.
- Docs and comments describe the method in general terms: no results from
  specific shots or cases.
- Plain, direct prose in docs and commit messages.
