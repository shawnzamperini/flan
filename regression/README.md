# Regression cases

Cases that run Flan on a simple plasma and check the result against an
analytic (fluid) result. They show that a piece of physics works, and
catch changes to the code that alter it.

The cases named `reg_*.cpp` are the current ones. Each is a normal Flan
input file, and everything about the case is set in it: the plasma comes
from the slab background (`bkg_source = "slab"`), a box of uniform plasma
in a uniform magnetic field with nothing hard-coded. The other input
files in this directory are older and are not part of this set.

| Case | Checks | Fluid result |
|---|---|---|
| `reg_friction` | Tungsten at rest in a flowing plasma is dragged up to the flow speed | Spitzer slowing-down time |
| `reg_thermal_force` | Tungsten in a plasma with parallel heat fluxes drifts against them | Drift where friction cancels 2.65 Z² ∇T_i + 0.71 Z² ∇T_e |

Each case has variants (other charge states, collision sub-cycling off,
the GPU, one heat flux at a time), chosen by a command-line argument and
listed at the top of its file.

## Building

Build and install Flan first (`make install` in its build directory),
then:

```
cmake -S <flan>/regression -B <build dir>
cmake --build <build dir>
```

If CMake finds a different MPI than Flan was built with, pass the one Flan
uses, e.g. `-DMPI_CXX_COMPILER=<path to mpicxx>`.

## Running

In the directory where the output should go:

```
for v in z10 z10_nosub z5 z10_gpu; do <build dir>/reg_friction $v; done
for v in z10 z10_ions z10_elec z10_nosub z10_off z5; do
    <build dir>/reg_thermal_force $v
done
```

Leave out `z10_gpu` without a GPU. A run takes from seconds to a few
minutes on a laptop. Each writes `<case>_<variant>.nc`.

## Checking

```
python <flan>/regression/check_collisions.py <output directory>
```

It reads every `reg_*` output file it finds there, takes all the inputs
from the files themselves, prints a table and writes `check_friction.png`
and `check_thermal_force.png`. The exit status is 1 if a check failed.

Each variant is checked against two things:

- **The fluid result, to 10%.** This is the physics check.
- **What the collision model should give at the run's time step, to 3%**
  plus the statistical error. Flan's collision model (Nanbu's) sits a few
  percent below the fluid result at any finite time step, by an amount
  that can be calculated (see the thermal force page of the docs). This
  tighter check catches changes the 10% check would let through.

Some variants are meant to miss the fluid result: those with collision
sub-cycling off, which show what sub-cycling corrects, and
`reg_thermal_force z10_elec`, where the default cap on sub-steps is too
small for the electrons (see the note at the top of
`reg_thermal_force.cpp`). The script recognizes them from the inputs:
most of the force being checked comes from a species whose sub-steps are
still more collisional than the target. Their table entry reads `miss
(expected)` and doesn't count as a failure; they are still held to the
model's expectation.

The numbers in the table (slowing-down times, drifts) can also be
compared between versions of the code.
