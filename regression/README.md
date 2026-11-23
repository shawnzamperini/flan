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
| `reg_prompt_redep` | Tungsten sputtered off a wall in a field parallel to it, ionizing to W+ and gyrating, comes back to the wall in its first gyration | Fussmann's prompt redeposition, ⟨r⟩(p) averaged over the exponential distribution of ionization lengths (`fussmann.py`) |
| `reg_sheath_angles` | D⁺ test ions fall from 10 ρ_s through the `"borodkina"` sheath field onto a wall, with the field at 45° or 60° to the normal | Energy conservation (exactly), and Borodkina and Tsvetkov's impact angle distributions (Fig. 2a of their 2015 paper) |

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
export OPENADAS_ROOT=<directory with scd50/ and acd50/>
for v in p01 p03 p1 p3 p10 b4 full drop0 brooks_ref; do
    <build dir>/reg_prompt_redep $v
done
for m in brooks borodkina; do for a in 2 5 10 20; do
    <build dir>/reg_prompt_redep ${m}_a$a
done; done
for v in a45 a60; do <build dir>/reg_sheath_angles $v; done
```

Leave out `z10_gpu` and `gradb_gpu` without a GPU. The collision runs
take from seconds to a few minutes each on a laptop, the drift runs about
a second. The prompt redeposition runs take about one (p01) to several
(p10) minutes each on 16 cores. Each writes `<case>_<variant>.nc`.
The sheath angle runs take about 2 (a45) and 4 (a60) minutes on 16
cores, and write `reg_sheath_angles_<variant>.txt`.

`reg_prompt_redep` reads the D on W yield table from `data/sputtering/` in
Flan's directory (set when it's built) and the OpenADAS rates from
`OPENADAS_ROOT`.

## Checking

```
python <flan>/regression/check_collisions.py <output directory>
python <flan>/regression/check_drifts.py <output directory>
python <flan>/regression/check_redeposition.py <output directory>
python <flan>/regression/check_sheath_angles.py <output directory>
```

Each reads the output files of its cases found there, takes all the
inputs from the files themselves, prints a table and writes plots
(`check_friction.png`, `check_thermal_force.png`, `check_drifts.png`,
`check_redeposition.png`, `check_sheath_angles.png`). The exit status is
1 if a check failed.

### Prompt redeposition

Fussmann's model gives the fraction of sputtered atoms that come back to
the wall in their first gyration after ionizing, for a field parallel to
the wall, cosine emission, and no electric field. It depends only on
p = e B / (m n_e S0), the mean ionization length over the gyroradius. For
an atom that ionizes after x gyroradii it is 1 / (1 + x²); Flan samples
ionization as a Poisson process, so the comparison is with that averaged
over the exponential distribution of x, ⟨r⟩(p) (`fussmann.py`, which
also checks 1 / (1 + x²) from the geometry by Monte Carlo).

The runs meet the model's assumptions with a few settings:

- The field is at 0.01 degrees to the wall, so the plasma reaches it, but
  an ion moves only about 1e-3 of a gyroradius along B toward the wall in
  a gyration.
- Ionization stops at W+ (`imp_iz_max_charge = 1`). Otherwise W+ would
  ionize again about 2π (S1 / S0) / p times per gyration, close to once at
  p = 1, which the model leaves out. The `full` variant allows it, to show
  how much it matters; it's printed but not checked.
- Collisions, the thermal force and the electric field are off.

The script works out p from the inputs saved in each file, with S0 from
the same OpenADAS file Flan used, interpolated the same way. A run passes
if Flan's prompt fraction (returns within one W+ gyroperiod of first
ionizing, over the source) is within 4 standard errors plus 0.01 of
⟨r⟩(p).

The cases with the sheath field on (`sheath_model`, `"off"` in the
others):

- `drop0` turns the field on with a drop of 1e-4 Te (`sheath_drop`).
  That is Fussmann's case again, and it is checked the same way.
- `<model>_a<angle>` put B at 2° to 20° to the wall with either profile.
  The script prints the prompt fraction and the redeposited ions' mean
  charge, impact energy and angle, and plots them against the angle
  (`check_redeposition_sheath.png`). They aren't checked; the gap between
  the profiles is the sensitivity to the profile's shape.
- `brooks_ref` is Brooks' reference case (D-T at 75 eV, 1e20 m⁻³,
  B = 5 T at 87° from the normal, f_D = 0.25 and a 3 Te drop set with
  `sheath_fd` and `sheath_drop`, collisions on). It is printed beside his
  numbers, not checked. Flan's W⁰ ionizes within about 0.06 mm (OpenADAS),
  several times nearer the wall than Brooks' 0.40 mm, so deep in the
  magnetic presheath: it gains less of the drop and comes back sooner, at
  a lower charge and nearer the normal. That difference comes from the
  ionization rates, not the sheath.

### Sheath impact angles

D⁺ test ions start 10 ρ_s from a slab wall with a Maxwellian at Ti
shifted by c_s along B (only those moving toward the wall), and fall
through the `"borodkina"` sheath field with no collisions, taking the
same steps as Flan's main loop, sub-steps near the wall included. Two
checks:

- **Energy conservation, to 1e-3 of the drop.** Each ion should hit with
  its starting energy plus e(φ(start) − φ_wall), whatever its orbit. This
  is the exact check of the field and the near-wall push.
- **The impact angle distribution** against Borodkina and Tsvetkov's
  model curves, read off their figure: the mean within 5° and the
  standard deviation within 3°. Their paper doesn't say how the
  Maxwellian is shifted, and the shift moves the mean by several degrees
  (29° to 39° at 60° for 0.5 to 2 c_s), so the mean's tolerance is wide.
  At 45° Flan's distribution lies on theirs; at 60° it has the same shape
  and width, about 4° nearer the normal.

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
