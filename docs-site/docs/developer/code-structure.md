# Code Structure

Flan is a C++ library, `libflan.so`. A simulation is a small C++ program, the
*input file*, that fills in an `Inputs` map of options and calls `flan()`.
Everything else happens inside the library. This page describes how the
library is organized, how a run flows through it, and how the CPU and GPU
versions of the code fit together. The [API Reference](../api/namespaces.md)
has the documentation for every function, generated from the comments in the
source code.

## Repository layout

| Directory | Contents |
|---|---|
| `src/` | The library source. One `.cpp` file per part of the code (options, background, Boris push, collisions, ...). |
| `src/cuda/` | The GPU (CUDA) versions of the main loop steps, plus the code that copies data to and from the GPU. Only built with CUDA enabled. |
| `include/` | Headers for `src/`. |
| `include/cuda/` | Headers for `src/cuda/`, and the small inline helpers used inside kernels (interpolation, indexing, ...). |
| `tests/` | Unit tests (GoogleTest), run with `ctest`. |
| `regression/` | Input files for the test cases (drifts, friction force) and a blank input file to start from. |
| `python/` | Helper scripts: reading Gkeyll data with postgkyl (`read_gkyl.py`), generating the Nanbu table (`calc_nanbu_A.py`) and plotting (`flan_plots.py`). |
| `scripts/` | Case setup (`flan_setup.sh`) and job submission scripts. |
| `docs-site/` | This website. |

`src/tmp/` and `include/tmp/` hold old code kept for reference. They aren't
built.

## How a run works

An input file looks like this (see `regression/blank_input_file.cpp`):

```cpp
int main(int argc, char* argv[])
{
	MPI_Init(&argc, &argv);
	Inputs inpts {create_inputs()};  // Fill in the options
	flan(inpts);
	MPI_Finalize();
}
```

MPI must be initialized before `flan()` is called, and every rank calls it.
`flan()` (`src/flan.cpp`) then runs the simulation in five steps:

1. **Load the options.** `Options::load_options` starts from the defaults and
   overwrites them with the values in `inpts`. An unrecognized option name,
   a value of the wrong type or an invalid value stops the run.
2. **Read the background.** `Background::read_bkg` loads the plasma on rank 0,
   from a Gkeyll simulation (`Gkyl::read_gkyl`) or one of the built-in test
   backgrounds (`Test::read_test`). It's then broadcast to every rank.
3. **Follow the impurities.** `ImpurityTransport::follow_impurities` follows
   this rank's share of the particles and records their statistics.
4. **Combine the statistics.** `Impurity::reduce_stats` adds up the
   statistics from every rank on rank 0, which then calculates the density,
   average velocities, average charge and average Nanbu `s` from them.
5. **Save the results.** `SaveResults::save_results` writes everything to
   `<case_name>.nc` on rank 0, and the timing summary is printed.

## The main loop

`follow_impurities` sets up the particle slots, the OpenADAS rate
coefficients and the random number generators, copies everything to the GPUs
if they're on, and then runs `main_loop` (`src/impurity_transport.cpp`).
Before the loop starts, the slots are filled, the starting cells are found
and the velocities are pushed back half a time step for the leapfrog scheme.
Then each pass through the loop runs these steps, in order:

| Step | Wrapper | CPU version | GPU version |
|---|---|---|---|
| Ionization/recombination (if on) | `ioniz_recomb_wrapper` | `ioniz_recomb_cpu` | `OpenADAS::ioniz_recomb_gpu` |
| Collisions (if on) | `collision_wrapper` | `collision_cpu` | `Collisions::collision_gpu` |
| Boris push | `boris_wrapper` | `Boris::update_velocity_cpu` | `Boris::update_velocity_gpu` |
| Step | `step_wrapper` | `step_cpu` | `step_gpu` |
| Bounds check | `check_bounds_wrapper` | `Boundary::check_bounds_cpu` | `Boundary::check_bounds_gpu` |
| Fill slots | `fill_slots_wrapper` | `Slots::fill_slots_cpu` | `Slots::fill_slots_gpu` |
| Find cell | `find_containing_cell_wrapper` | `find_containing_cell_cpu` | `find_containing_cell_gpu` |
| Record statistics | `record_stats_wrapper` | `Impurity::record_stats_cpu` | `ImpurityStats::record_stats_gpu` |

The loop ends once every slot is dead and no particles are left to launch.
The physics of each step is described in
[Impurity Transport](../physics/impurity-transport.md),
[Collisions](../physics/collisions.md) and
[Ionization/Recombination](../physics/ionization-recombination.md).

### Particle slots

The particles being followed are held in a fixed number of *slots*
(`Slots::Slots`, `src/slots.cpp`). There are `min(slot_cap, imp_num)` of
them. When a particle dies (it's absorbed at a boundary, or exceeds
`imp_max_time`), its slot is refilled with a new particle in the fill slots
step, until all of this rank's particles have been launched. The slots are
stored as a *structure of arrays*: one array per quantity (`t`, `x`, `vX`,
`q`, `state`, ...), indexed by slot. This is what lets each GPU thread work on
one slot with fast, contiguous memory access.

The particles are split as evenly as possible between the MPI ranks (the last
rank gets the remainder), and then between the GPUs on a rank (GPU 0 gets the
remainder).

## CPU and GPU code

Every main loop step has a CPU version and a GPU version, and a *wrapper* that
picks between them based on the `use_gpu` option. The wrappers are in
`src/impurity_transport.cpp`, and they all look like this:

```cpp
void step_wrapper(...)
{
#ifdef USE_CUDA
	if (opts.use_gpu_int() > 0)
	{
		// Defined in impurity_transport.cu
		step_gpu(slots_d, opts.imp_time_step(), opts_d, warnings_d.max_time);
		return;
	}
#endif

	// Defined above
	step_cpu(slots, bkg, opts, opts.imp_time_step(), max_time_ctr);
}
```

The physics is the same either way. The two versions are kept as close to
each other as possible, so a change to one usually needs the same change in
the other.

### CPU (OpenMP)

The CPU versions loop over the slots with `#pragma omp parallel for`. Each
OpenMP thread has its own random number generator.

### GPU (CUDA)

The GPU versions are split into a *kernel* (e.g., `step_kernel`) that runs one
thread per slot, and a host function (e.g., `step_gpu`) that launches it in
blocks of 256 threads. Launches are asynchronous, except in `DEBUG` builds,
where each one waits for the kernel to finish and checks for errors. Each
slot has its own random number generator (`init_rngs.cu`).

Kernels can't use the host-side classes, so the data they need is copied to
the GPU into plain structs of device pointers:

| Host side | Device side | Copied with |
|---|---|---|
| `Background::Background` | `Background::BackgroundDevice` | `Background::to_device` |
| `Slots::Slots` | `Slots::SlotsDevice` | `Slots::to_device` |
| `Impurity::Statistics` | `ImpurityStats::StatisticsDevice` | `Statistics::to_device` |
| `OpenADAS::OpenADAS` | `OpenADAS::OpenADASDevice` | `OpenADAS::to_device` |
| `Options::Options` | `Options::OptionsDevice*` (only the options kernels use) | `Options::to_device` |

Each has a matching `free_*` function. The background's times, cell centers,
cell edges and bounds are small and read by nearly every kernel, so they go in
GPU *constant memory* instead (`src/cuda/device_constants.cu`). Constant memory
is shared by everything on a GPU, so only one background can be on each GPU.

With several GPUs, `follow_impurities` copies everything to each GPU and runs
`main_loop` on each in its own `std::thread` (`gpu_worker`). Counters that the
GPUs add to, like the ionization/recombination and `imp_max_time` warnings,
are kept on each GPU and summed on the host after the threads finish.

### Fields at a particle

Interpolating the background at a particle's position and time is the most
common operation in the main loop. On the GPU it's done by building a
*stencil* once per particle (`Interpolate::build_stencil_4d`, in
`include/cuda/interpolate.cuh`) and reusing it for every field. The fields
are interpolated trilinearly in space and linearly in time.

## Data layout

The background and statistics arrays are `Vectors::Vector4D`s with dimensions
(t, x, y, z), and the geometry arrays are `Vectors::Vector3D`s with dimensions
(x, y, z). Both store the data as one flattened vector in row-major order, so
the last index varies fastest. The same layout is used on the GPU
(`Utilities::calc_4d_index_cuda`) and in the netCDF output, so the data is
copied as is between them.

## Namespaces

| Namespace | Contents |
|---|---|
| `Options` | The `Options` class holding every input option, and loading them from the input file. |
| `Background` | The `Background` class holding the plasma and geometry, and `read_bkg`. |
| `Gkyl` | Reading a Gkeyll background. |
| `Test` | The built-in test backgrounds. |
| `Slots` | The particle slots, and creating new particles. |
| `ImpurityTransport` | `follow_impurities`, the main loop and the CPU/GPU wrappers. |
| `Boris` | The Boris push. |
| `Collisions` | The Nanbu collision model. |
| `OpenADAS` | Reading the OpenADAS rate coefficients, and ionization/recombination. |
| `Boundary` | The boundary conditions. |
| `Impurity` | The `Statistics` class holding the results, recording them on the CPU, and `reduce_stats`. |
| `ImpurityStats` | Recording the statistics on the GPU, and the device-side `StatisticsDevice`. |
| `SaveResults` | Writing the netCDF output file. |
| `Interpolate`, `Indices`, `Utilities` | Interpolation, cell lookup and other helpers. |
| `Vectors` | `Vector3D` and `Vector4D`. |
| `Errors` | `fatal_error` and `check_cuda` (see below). |
| `Timer`, `Constants`, `Nanbu` | Timing, physical constants and the tabulated Nanbu A(s). |

## Errors

An error that would give wrong results if the run carried on stops it with
`Errors::fatal_error(msg)` (`src/errors.cpp`). It prints the message with the
MPI rank and calls `MPI_Abort`, so every rank stops, not just the one that hit
the error. `Errors::check_cuda` does the same for a failed CUDA call. Device
code can't stop the run itself. A kernel that finds a problem the run can
continue past counts it instead, and the host reports the count at the end.

## Building and testing

Flan is built with CMake. Every `src/*.cpp` file is compiled into the library,
and `src/cuda/*.cu` too when CUDA is enabled (`-DENABLE_CUDA=ON`), which also
defines `USE_CUDA`. CPU-only builds get placeholder versions of the GPU
functions in the `.cpp` files, inside `#ifndef USE_CUDA`, that stop the run if
they're ever called.

The unit tests in `tests/` are built with the library. The `*_cpu` tests are
always built, and the `*_cuda` tests only with CUDA. Run them with `ctest`
from the build directory. The regression inputs in `regression/` are full
simulations of the test backgrounds, compiled against the library like any
other input file.

## Adding an input option

An option has to be added in several places, which don't check each other:

1. `include/options.h` and `src/options.cpp`: the member with its default, a
   getter and a setter. A string option also needs its control integer, set in
   the setter, and the setter called in `initialize_control_ints`.
2. `src/load_options.cpp`: the option name in `load_input`, so it can be set
   from an input file.
3. `src/save_results.cpp`: the option in `save_input_options`, so it's saved
   in the output file.
4. If a kernel needs it: a field in `OptionsDevice` (`include/options_device.h`),
   set in `Options::to_device` (`src/cuda/options.cu`). A field that isn't set
   there is zero on the GPU.
5. `docs-site/docs/input-options.md`: the option, its type, default and
   description.
