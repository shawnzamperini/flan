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
| `reg_sheath_angles` | D⁺ test ions fall from 10 ρ_s through the `"borodkina"` sheath field onto a wall, with the field at 45° or 60° to the normal, pushed by the main loop's routines or (`_sampler`) by `Sheath::sample_impacts` | Energy conservation (exactly), Borodkina and Tsvetkov's impact angle distributions (Fig. 2a of their 2015 paper), and the sampler against the main loop |
| `reg_self_sputter` | Tungsten landing on the wall sputters more tungsten, with a constant yield Y0 (0.5, 1.25) | The weight in each generation, W_{g+1} = Y0 L_g, and the total S f₀ Y0 / (1 − f₁ Y0), with the self-sputtered atoms coming back as often as the D-sputtered ones |

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
for v in a45 a60 a45_sampler a60_sampler; do <build dir>/reg_sheath_angles $v; done
python <flan>/regression/make_self_tables.py
for v in y05 y125 w_te20 w_te50; do <build dir>/reg_self_sputter $v; done
```

Leave out `z10_gpu` and `gradb_gpu` without a GPU. The collision runs
take from seconds to a few minutes each on a laptop, the drift runs about
a second. The prompt redeposition runs take about one (p01) to several
(p10) minutes each on 16 cores. Each writes `<case>_<variant>.nc`.
The sheath angle runs take about 2 (a45) and 4 (a60) minutes on 16
cores, and a few seconds for the `_sampler` variants, and write
`reg_sheath_angles_<variant>.txt`.

`reg_prompt_redep` reads the D on W yield table from `data/sputtering/` in
Flan's directory (set when it's built) and the OpenADAS rates from
`OPENADAS_ROOT`.

`reg_self_sputter` needs the same, plus the constant-yield tables
`self_y05.h5` and `self_y125.h5` in the working directory, which
`make_self_tables.py` writes, and `data/sputtering/W_on_W.h5` for the
`w_te<Te>` variants. y05 takes about 4 minutes on 16 cores.

## Checking

```
python <flan>/regression/check_collisions.py <output directory>
python <flan>/regression/check_drifts.py <output directory>
python <flan>/regression/check_redeposition.py <output directory>
python <flan>/regression/check_sheath_angles.py <output directory>
python <flan>/regression/check_self_sputter.py <output directory>
```

Each reads the output files of its cases found there, takes all the
inputs from the files themselves, prints a table and writes plots
(`check_friction.png`, `check_thermal_force.png`, `check_drifts.png`,
`check_redeposition.png`, `check_sheath_angles.png`,
`check_self_sputter.png`). The exit status is
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

### Self-sputtering

With `sput_self = "on"`, a W atom or ion landing on the wall sputters Y
atoms and the particle carries on as one of them with its weight times Y.
The output tallies the weight carried into each generation, W_g (W₀ = S,
the D⁺ source), and the weight landing on the wall from it, L_g, so
f_g = L_g / W_g is generation g's return fraction.

`y05` and `y125` are `reg_prompt_redep`'s p1 case (W⁺ only, B at 0.01° to
the wall, no electric field) with a constant yield, 0.5 or 1.25, from
`make_self_tables.py`. With a constant yield every particle of a
generation has the same weight (S / imp_num times Y0^g, or the roulette
weight), so the script can count them, and checks:

- W_{g+1} = Y0 L_g: exactly where no roulette is played (and the counts
  are whole numbers), and within 4 standard deviations of the binomial
  roulette where it is.
- The self-sputtered total, its sum over the elements and the redeposited
  total agree with the generation tallies.
- f₁ (generations 1 and up) equals f₀ within 4 standard errors plus 0.02.
  Prompt redeposition doesn't depend on an atom's speed, so atoms
  sputtered by W, which are faster, should come back as often as those
  sputtered by D⁺.
- The self-sputtered rate is S f₀ Y0 / (1 − f₁ Y0), and the
  multiplication T / (S + T) is f₁ Y0, to 2%.

`w_te<Te>` use the RustBCA W on W table at Te = Ti = 20 or 50 eV with
full ionization and the sheath on, and are printed, not checked: the
self-sputtered rate, the multiplication, the landing ions' mean charge,
energy and yield, and the yield at the energies the charges fall through.

