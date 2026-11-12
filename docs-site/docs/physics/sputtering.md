# Sputtering

Ions from the background plasma that hit a wall knock atoms out of it. This
is physical sputtering, and it's the main source of impurities such as
tungsten from divertor targets. How many atoms come out per incident ion is
the sputtering yield $Y$, which depends on the projectile, the wall
material, and the ion's energy and angle of incidence.

!!! note "Work in progress"
    Flan can read and interpolate sputtering yield tables, which is what
    this page describes. Launching impurities from the walls of a SOLEDGE3X
    background using these yields is being added next, and this page will
    be extended to cover it. Until then, the yield tables aren't used in a
    run.

## Where the yields come from

The yields are computed with [RustBCA](https://github.com/lcpp-org/RustBCA)
(Drobny and Curreli, 2021), a binary collision approximation (BCA) code.
A BCA code follows the incident ion and the target atoms it sets moving
through a series of two-body collisions, and counts the target atoms that
leave the surface with enough energy to escape.

RustBCA is installed into the `flan` conda environment when the
environment is set up (see [Installation](../installation/installation.md#rustbca)),
pinned to version 3.0.0 (commit `51bee88`). The tables are made ahead of
time with `python/rustbca_yields.py`, which calls RustBCA's
`sputtering_yield` function, and are kept in the repository in
`data/sputtering/`. Running Flan doesn't need RustBCA, only making a new
table does.

`sputtering_yield` simulates a single-element, semi-infinite target with a
flat surface. It uses RustBCA's default models, which can't be changed
through that function:

| Setting | Model |
|---|---|
| Interatomic potential | Kr-C |
| Electronic stopping | Lindhard-Scharff, nonlocal (`LOW_ENERGY_NONLOCAL`) |
| Mean free path | Liquid (amorphous target) |
| Scattering integral | Mendenhall-Weller quadrature |
| Surface binding | Average (`AVERAGE`), with the target's surface binding energy $E_s$ |
| Bulk binding | Individual, with bulk binding energy $E_b$ |

These are written into each table file along with the RustBCA version, so
it's always known how a table was made.

### Angle convention

The angle of incidence is measured from the surface normal, in degrees: 0
is head-on and 90 is grazing.

### Statistics

The yield at each energy and angle is the number of sputtered atoms
divided by the number of ions run, $N$. Its error is the counting error,
$\sqrt{Y N}/N$. Near threshold, $Y$ can be $10^{-4}$ or smaller, so a fixed
$N$ is either too noisy there or very slow everywhere else. Instead, each
point is run twice:

1. A pilot run with $N_\text{min}$ ions ($10^4$ by default) estimates $Y$.
2. A run with enough ions to count about 2000 sputtered atoms (about 2%
   error), $N = 2000 / Y$, but no fewer than $N_\text{min}$ and no more than
   $N_\text{max}$ ($2\times10^6$ by default). If the pilot run sputtered
   nothing, $N_\text{max}$ ions are used.

RustBCA runs every ion with the same seed, and ion $i$ always uses random
stream $i$, so running more ions repeats the first ones rather than adding
to them. This is why the second run replaces the pilot run rather than
adding to it. The seed is set with `--seed` (default 0), so a table can be
remade exactly.

The energies are run from highest to lowest. Once an angle sputters nothing
with $N_\text{max}$ ions, that energy is below threshold at that angle, and
so are the lower energies. They're set to $Y = 0$ without running them
(their `num_samples` is 0). This saves running millions of ions just to get
zero.

### Energies near threshold

Just above threshold the yield rises very steeply. For D on W it goes up
about fivefold between neighboring energies of a 41-point grid from 10 eV
to 10 keV (19% apart), and interpolating that grid was off from RustBCA by
17% at 280 eV and by a factor of 3.5 at 230 eV. Interpolating $\log Y$
instead didn't help. This is the range that matters most for D on W, since
D ions hit a divertor target with roughly $3T_e + 2T_i$, a few hundred eV
at typical temperatures.

So after the main grid is run, energies about 5% apart are added across
the threshold region: from the highest energy where every angle sputtered
nothing up to 3 times the lowest energy where any angle sputtered
something (`--fine-step` and `--fine-factor`, or `--no-refine` to turn
this off). For D on W this adds 25 energies between 197 and 672 eV.

### Remaking or extending a table

Every point is independent, and RustBCA is seeded the same way for each,
so a point comes out the same in any table made with the same materials,
seed, `--target-counts`, `--n-min` and `--n-max`. `--reuse FILE` copies the
points of an existing table instead of running them again, e.g., to add
energies or angles to a table without remaking the rest. The script checks
that the existing table was made with the same settings, which are saved
as attributes of each table, and stops if not.

## Table files

Each projectile and target pair has its own HDF5 file,
`data/sputtering/{projectile}_on_{target}.h5`, e.g., `D_on_W.h5` for
deuterium on tungsten.

| Dataset | Shape | Units | Contents |
|---|---|---|---|
| `/energy` | (nE) | eV | Incident energies, log spaced |
| `/angle` | (nA) | degrees | Incident angles from the surface normal |
| `/yield` | (nE, nA) | atoms/ion | Sputtering yield |
| `/yield_err` | (nE, nA) | atoms/ion | Counting error of the yield |
| `/num_samples` | (nE, nA) | | Number of ions run for each point |
| `/projectile/Z`, `/projectile/mass_amu` | scalar | , amu | Projectile atomic number and mass |
| `/projectile/Es`, `/projectile/Ec` | scalar | eV | Projectile surface binding and cutoff energies |
| `/target/Z`, `/target/mass_amu` | scalar | , amu | Target atomic number and mass |
| `/target/n` | scalar | m$^{-3}$ | Target number density |
| `/target/Es`, `/target/Eb`, `/target/Ec` | scalar | eV | Target surface binding, bulk binding and cutoff energies |

The root group's attributes record where the table came from: the
projectile and target names, the RustBCA version and commit, the date, the
command that made it, the RustBCA models above, and the sampling settings
(`seed`, `target_counts`, `n_min`, `n_max`, `refine`, `fine_step` and
`fine_factor`). Flan only
reads datasets, so everything it needs is a dataset rather than an
attribute.

A table can be looked at with `h5dump` or h5py, e.g.,

```python
import h5py
with h5py.File("data/sputtering/D_on_W.h5") as f:
    print(dict(f.attrs))
    print(f["energy"][:], f["angle"][:], f["yield"][:])
```

### Tables in the repository

| File | Projectile | Target | Energies | Angles |
|---|---|---|---|---|
| `D_on_W.h5` | D | W | 66: 41 from 10 eV to 10 keV, plus 25 from 197 to 672 eV | 20, 0 to 85° every 5°, 87° and 89° |

The D on W yields, as a function of energy at a few angles and as a
function of angle at a few energies (made with `--plot-only` and
`--save-plots`). The threshold is about 200 eV at every angle. Above a few
keV the yield rises toward grazing incidence, peaking around 85°, and
near threshold it's highest at moderate angles.

![D on W sputtering yield vs. energy](../assets/images/sputtering/D_on_W_energy.png)

![D on W sputtering yield vs. angle](../assets/images/sputtering/D_on_W_angle.png)

Interpolating this table agrees with RustBCA run at energies and angles
between the table's points to within the statistical error: within about
5% above 280 eV, and within the (large) counting error of RustBCA itself
just above threshold, where the yield is below $10^{-5}$.

## How Flan uses a table

A table is read by `Sputtering::YieldTable` (`src/sputtering.cpp`), which
checks it as it's read: the run stops with a message if a dataset is
missing, the yield has the wrong shape, there are fewer than two energies
or angles, the energies aren't positive and increasing, the angles aren't
increasing within 0 to 90°, or a yield is negative.

The yield at energy $E$ and angle $\theta$ is bilinearly interpolated,
linear in $\log_{10} E$ and in $\theta$. $Y$ itself is interpolated rather
than $\log Y$, because $Y = 0$ below threshold. Outside the table:

- Below the lowest energy, $Y = 0$. Tables start below threshold, so this
  is the right answer, not an extrapolation.
- Above the highest energy, the yield at the highest energy is used.
- Angles are clamped to the table's range, so 90° uses the largest angle
  in the table (89° for the default angles).

## Making a new table

Activate the `flan` environment, which has RustBCA, and run the script
with the projectile and target:

```bash
conda activate flan
python python/rustbca_yields.py D W
```

This writes `data/sputtering/D_on_W.h5`. It prints each energy as it
finishes. The `D_on_W.h5` table took about 3 hours on 20 cores (RustBCA
runs the ions in parallel on every core): 70 minutes for the main grid and
100 minutes for the energies near threshold, where most points need the
full $N_\text{max}$ ions. Before a full run, try
`--quick`, which uses a coarse grid and few ions and takes seconds:

```bash
python python/rustbca_yields.py D W --quick --out test.h5 --plot
```

The options are:

| Option | Default | Description |
|---|---|---|
| `projectile`, `target` | `D`, `W` | Names from the script's `MATERIALS` table. A target needs a density and bulk binding energy. |
| `--out` | `data/sputtering/{projectile}_on_{target}.h5` | Output file |
| `--emin`, `--emax` | 10, 10000 | Energy range (eV) |
| `--ne` | 41 | Number of log-spaced energies |
| `--angles` | 0 to 85 every 5, 87, 89 | Angles from the normal (degrees), between 0 and 90 |
| `--target-counts` | 2000 | Sputtered atoms to aim for at each point |
| `--n-min` | 10000 | Ions in each pilot run, and the fewest per point |
| `--n-max` | 2000000 | Most ions per point |
| `--seed` | 0 | RustBCA seed |
| `--no-refine` | | Don't add energies near threshold |
| `--fine-step` | 0.05 | Fractional spacing of the energies added near threshold |
| `--fine-factor` | 3 | Energies are added up to this times the lowest energy that sputters anything |
| `--reuse FILE` | | Copy points from an existing table made with the same settings instead of running them |
| `--target-Es` | | Override the target's surface binding energy (eV) |
| `--quick` | | Coarse grid (9 energies, 4 angles) and few ions, for testing |
| `--plot` | | Plot the table after making it |
| `--plot-only FILE` | | Only plot an existing table |
| `--save-plots PREFIX` | | Save the plots as `PREFIX_energy.png` and `PREFIX_angle.png` instead of showing them |

### Adding a projectile or target

The materials are defined in the `MATERIALS` dictionary at the top of
`python/rustbca_yields.py`. Add an entry with:

| Key | Units | Needed for | Description |
|---|---|---|---|
| `Z` | | both | Atomic number |
| `m` | amu | both | Mass |
| `Es` | eV | both | Surface binding energy. For a target, usually the heat of sublimation. |
| `Ec` | eV | both | Cutoff energy, below which RustBCA stops following a particle |
| `Eb` | eV | target | Bulk binding energy (often 0 for metals) |
| `n` | m$^{-3}$ | target | Number density |

RustBCA's `scripts/materials.py` (in `$FLANSOFT/RustBCA`) has values for
many elements, but calls them example values, so check them. The values in
use now are:

| Material | Z | m (amu) | $E_s$ (eV) | $E_c$ (eV) | $E_b$ (eV) | n (m$^{-3}$) |
|---|---|---|---|---|---|---|
| D | 1 | 2.014 | 1.5 | 0.1 | | |
| W | 74 | 183.84 | 8.79 | 3.0 | 0 | $6.306\times10^{28}$ |

$E_s$ for W is its heat of sublimation, the same binding energy that will
be used for the energy distribution of sputtered atoms (Thompson, 1968). RustBCA's `materials.py` uses
8.9 eV. Near threshold the yield is sensitive to $E_s$.

Then make the table, look at the plots, check it against reference data
(see [Validation](#validation)), add it to the
[table of files above](#tables-in-the-repository), and commit it with the
command that made it in the commit message.

## Validation

!!! note "To do"
    Compare the D on W yields with Eckstein's tabulated values (Eckstein,
    2007; Eckstein, 2002) and the Bohdansky-Eckstein fit, at normal
    incidence and as a function of angle, and add the plots here. Expected
    features: a threshold of roughly 200 eV at normal incidence, and at
    high energy a yield that rises toward grazing incidence before falling
    near 90°.

<!--
Plots to add, made with `python python/rustbca_yields.py --plot-only
data/sputtering/D_on_W.h5 --save-plots ...` plus the reference data:

![D on W yield vs. energy, RustBCA and Eckstein](../assets/images/sputtering/D_on_W_energy_eckstein.png)

![D on W yield vs. angle, RustBCA and Eckstein](../assets/images/sputtering/D_on_W_angle_eckstein.png)
-->

## Limitations

- **Pure, fresh targets.** Each table is for a single-element target with
  nothing implanted in it. A real wall holds deuterium and other
  impurities near its surface, and can be rough, both of which change the
  yield.
- **Physical sputtering only.** Chemical sputtering, sputtering of loosely
  bound surface atoms (ad-atoms) and evaporation aren't included.
- **Threshold region.** Near threshold, the yield depends strongly on the
  surface binding model and $E_s$, and BCA codes are least accurate there.
  For D on W this is the energy range that matters most at typical
  divertor temperatures, so the uncertainty in the yield is largest where
  it counts.
- **No sputtering by the impurities themselves** (self-sputtering) or by
  other impurities in the plasma yet.

## References

- J. Drobny and D. Curreli, "RustBCA: A high-performance
  binary-collision-approximation code for ion-material interactions",
  *J. Open Source Softw.* **6**, 3298 (2021).
- W. Eckstein, "Sputtering yields", in *Sputtering by Particle Bombardment*,
  Topics Appl. Phys. **110**, 33 (2007).
- W. Eckstein, "Calculated sputtering, reflection and range values",
  IPP Report 9/132 (2002).
- M. W. Thompson, "The energy spectrum of ejected atoms during the high
  energy sputtering of gold", *Philos. Mag.* **18**, 377 (1968).
