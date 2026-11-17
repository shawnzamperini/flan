# Regression cases

Cases that run Flan on a simple plasma and check the result against an
analytic (fluid) result. They show that a piece of physics works, and
catch changes to the code that alter it.

The cases named `reg_*.cpp` are the current ones. Each is a normal Flan
input file, and everything about the case is set in it: the plasma comes
from the slab background (`bkg_source = "slab"`), a box of uniform plasma
in a simple magnetic field with nothing hard-coded. The other input
files in this directory are older and are not part of this set.

| Case | Checks | Fluid result |
|---|---|---|
| `reg_friction` | Tungsten at rest (or carbon at the ion temperature) in a flowing plasma is dragged up to the flow speed | Spitzer slowing-down time (between two Maxwellians for carbon) |
| `reg_thermal_force` | Tungsten (or carbon) in a plasma with parallel heat fluxes drifts against them | Drift where friction cancels Z² (β_i ∇T_i + 0.71 ∇T_e), with Stangeby's mass-dependent β_i (2.61 for W in D, 2.05 for C in D) |
| `reg_drifts` | Gyration, and the ExB, polarization, grad-B and curvature drifts of a single particle, without collisions | Cyclotron frequency and the guiding-center drift formulas |

Each case has variants (other charge states, collision sub-cycling off,
the GPU, one heat flux or one drift at a time), chosen by a command-line
argument and listed at the top of its file.

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
for v in z10 z10_nosub z5 z10_gpu c4; do <build dir>/reg_friction $v; done
for v in z10 z10_ions z10_elec z10_nosub z10_off z5 c4_ions; do
    <build dir>/reg_thermal_force $v
done
for v in gyration exb polarization gradb curvature gradb_gpu; do
    <build dir>/reg_drifts $v
done
```

Leave out `z10_gpu` and `gradb_gpu` without a GPU. The collision runs
take from seconds to a few minutes each on a laptop, the drift runs about
a second. Each writes `<case>_<variant>.nc`.

## Checking

```
python <flan>/regression/check_collisions.py <output directory>
python <flan>/regression/check_drifts.py <output directory>
```

Each reads the output files of its cases found there, takes all the
inputs from the files themselves, prints a table and writes plots
(`check_friction.png`, `check_thermal_force.png`, `check_drifts.png`).
The exit status is 1 if a check failed.

### Drifts

Every particle of a drift run starts the same way and there are no
collisions, so the mean velocity in the output is that of one particle: a
steady drift plus the gyration around it. The script separates the two
with a fit and compares the drift with theory to 2% of its size (0.5% for
the cyclotron frequency and speed in the gyration case).

### Collisions

Each variant is checked against two things:

- **The fluid result, to 10%.** This is the physics check.
- **What the collision model should give at the run's time step, to 3%**
  plus the statistical error. Flan's collision model (Nanbu's) sits a few
  percent below the fluid result at any finite time step, by an amount
  that can be calculated (see the thermal force page of the docs). This
  tighter check catches changes the 10% check would let through.

The model's expectation is worked out for an impurity that is slow next
to the ions, so the carbon variants (`c4`, `c4_ions`), which start at the
ion temperature, are checked against fluid only. For them the fluid
check is not circular: the heat-flux skew is fitted to the heavy limit,
and how the thermal force changes with the impurity's mass comes only
from the collisions.

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

### Effective thermal-force coefficients

After the checks, the script prints the thermal-force coefficients of the
`reg_thermal_force` runs: the β_i and α_e of F = Z² (β_i ∇T_i + α_e ∇T_e)
that Flan's drift corresponds to. At force balance the friction m_z ν_s
v_drift cancels the thermal force, so each coefficient is −m_z ν_s v_drift
/ (Z² ∇T). `z10_ions` gives β_i, `z10_elec` gives α_e, and the runs with
both heat fluxes give their sum. These numbers are for comparison only and
are not checked.

The friction rate ν_s is taken two ways. The fluid (Spitzer) rate gives the
coefficient that a fluid code would infer from Flan's drift. Flan's own
rate, measured by the `reg_friction` run with the same plasma, time step
and sub-cycling, isolates the force Flan actually applies. The two differ
by Nanbu's friction shortfall at a finite step. Because the friction and
the thermal force fall short by about the same factor, the drift is closer
to fluid than either force is on its own.

The table also lists three references:

- What the collision model should give at the run's step.
- The heavy-limit fluid value (2.65 or 0.71) that the heat-flux skew is
  fitted to. Getting it back is a consistency check, not an independent
  one.
- Stangeby's mass-dependent β_i(μ), with μ = m_z / (m_z + m_i).

∇T is not an input of the slab. It is the gradient that would carry the
run's heat flux in a collisional plasma (Braginskii's conductivities).
