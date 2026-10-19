# Running Flan

A Flan simulation is a *case*: a directory with an input file, which is a
small C++ program that sets the options and calls Flan. Running one takes
four steps:

1. Set up a case directory with `flan_setup.sh`.
2. Set the options in the input file.
3. Compile it with `make`.
4. Run it, from the command line or with a job script.

This page assumes Flan is installed (see [Installation](../installation/installation.md)).
Activate the `flan` conda environment before any of these steps:

```bash
conda activate flan
```

The environment has the paths Flan needs to read Gkeyll data. A Gkeyll run
without it stops with an error saying `READ_GKYL_PY` is not set.

## Set up a case

Create a case directory with the setup script, which is made when Flan is
built:

```bash
/path/to/flan/scripts/flan_setup.sh my_case
cd my_case
```

This creates `my_case/` with:

- `my_case.cpp`: the input file, a copy of `regression/blank_input_file.cpp`
  to start from.
- `CMakeLists.txt` and the `Makefile` made from it, which compile the input
  file against the installed Flan library.

The executable is named after the case, and it's built from the `.cpp` file
with the same name, so keep the two names the same.

## The input file

The input file sets the options in `create_inputs` and passes them to Flan:

```cpp
#include "flan.h"
#include "mpi.h"

Inputs create_inputs()
{
	Inputs inpts;

	inpts["case_name"] = "my_case";
	inpts["imp_num"] = 100000;
	inpts["imp_time_step"] = 1e-8;
	inpts["imp_atom_num"] = 74;
	inpts["imp_mass_amu"] = 183.84;
	// ... any other options

	return inpts;
}

int main(int argc, char* argv[])
{
	MPI_Init(&argc, &argv);
	Inputs inpts {create_inputs()};
	flan(inpts);
	MPI_Finalize();
}
```

Every option is listed on the [Input Options](../input-options.md) page, and
any option you don't set uses its default. The value has to have the
option's type: `inpts["imp_num"] = 1e5` is a double, so it stops the run for
an int option, while `inpts["imp_num"] = 100000` works. An option name that
isn't recognized, e.g., from a typo, also stops the run, before the
simulation starts.

The blank input file turns on ionization/recombination, which needs the
OpenADAS rate coefficient files. These aren't installed with Flan. Download
them from [open.adas.ac.uk](https://open.adas.ac.uk) and set `openadas_root`
to the directory they're in (see
[Ionization/Recombination](../physics/ionization-recombination.md) for the
layout Flan expects).

## Compile

```bash
make
```

This builds the executable `my_case`. Run `make` again after every change to
the input file, and after reinstalling Flan.

## Run from the command line

```bash
./my_case
```

This runs on one MPI rank, using as many OpenMP threads as OpenMP chooses,
usually one per core. To set the number of threads, e.g., to 8:

```bash
OMP_NUM_THREADS=8 ./my_case
```

To run on several MPI ranks, e.g., 2, with 8 threads each:

```bash
OMP_NUM_THREADS=8 mpiexec -n 2 ./my_case
```

Your system may use `mpirun` or `srun` instead of `mpiexec`. On a laptop or
workstation, one rank with all the cores as threads is usually best. Several
ranks are for spreading a run across the nodes of a cluster. Every rank
holds its own copy of the background plasma, so more ranks on one machine
use more memory, while more threads don't.

### GPUs

With a GPU build of Flan installed, set `inpts["use_gpu"] = "cuda"` to follow
the particles on the GPUs. Each rank uses every GPU it can see, so run one
rank per node (see [MPI Parallelism](../developer/mpi-parallelism.md)).

## What a run does

Flan prints a short banner with the version and the number of ranks and
threads, then loads the options and the background. While particles are
followed it prints progress updates (`print_interval` controls how many), and
at the end, any warnings and a summary of where the time went.

The results are saved to `<case_name>.nc`, a netCDF file with four groups:

| Group | Contents |
|---|---|
| `input` | Every input option used in the run. |
| `geometry` | The times, the grid and the geometry of the background. |
| `background` | The background plasma: density, temperatures, fields and flow. |
| `output` | The impurity results: counts, density `nz`, average velocities, average charge and Nanbu `s`. |

The 4D arrays have dimensions (time, x, y, z). The file can be read with any
netCDF reader, e.g., `netCDF4` in Python. The `flan` environment also has
`flan_plots`, with plotting helpers for Flan output:

```python
import flan_plots

fp = flan_plots.FlanPlots("my_case.nc")
fp.plot_frame_xy("nz", frame=0, z0=0.0)  # Impurity density at z = 0
```

### Gkeyll background files

The first time a case runs with a Gkeyll background, Flan reads the Gkeyll
data with postgkyl and saves it as `bkg_from_pgkyl_*` files in the case
directory. Later runs load those files directly, which is much faster.

!!! note "Delete the background files after changing the background"
    The saved files are reused no matter what the Gkeyll options say. If you
    change `gkyl_dir`, `gkyl_casename`, the frames or the species, delete the
    `bkg_from_pgkyl_*` files so they're made again. Flan stops if the number
    of frames doesn't match, but it can't catch every change.

## Example: the ExB drift test

The test cases in `regression/` are complete input files. To run the ExB
drift test, which follows 100 particles through constant electric and
magnetic fields:

```bash
/path/to/flan/scripts/flan_setup.sh test_exb
cd test_exb
cp /path/to/flan/regression/test_exb.cpp .
make
./test_exb
```

The copy replaces the blank input file, since it has the same name. The
results are in `exb.nc`, named after its `case_name`.

## Running with a job script

On a cluster, compile on a login node as above, then submit the run as a job.
There are example job scripts for:

- [Perlmutter](running-perlmutter.md)
- [Omega](running-omega.md)

## Troubleshooting

If a run stops with an `Error!` message, the message says what's wrong, e.g.,
an option name or value, a missing file or a background that doesn't match
the options. For anything else, contact Shawn Zamperini
([zamperinis@fusion.gat.com](mailto:zamperinis@fusion.gat.com)).
