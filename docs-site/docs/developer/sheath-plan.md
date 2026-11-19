# Sheath and Self-Sputtering Plan

This page is the development plan for adding the sheath in front of
material surfaces to Flan, so that prompt redeposition of sputtered
impurities is modelled properly, and then self-sputtering, where impurity
ions that hit the wall sputter more atoms. It builds on the sputtered
source (`imp_source_opt = "sputter"`, see [Sputtering](../physics/sputtering.md)),
which works on SOLEDGE3X walls and HERMES targets.

## Why

A sputtered W atom leaves the wall with a few eV, ionizes a fraction of a
millimeter from it, and starts to gyrate. Whether it comes back to the
wall in its first gyration (prompt redeposition) or gets away sets how
much W leaves the divertor. Two things decide it that Flan doesn't have
yet:

- **The sheath.** Within a few Debye lengths and a few ion sound
  gyroradii of the wall, a strong electric field pulls positive ions onto
  it. A W ion born inside that region is accelerated back to the wall
  before it can gyrate away, so prompt redeposition is higher than without
  it. Flan's background has no sheath: the fluid codes treat it as a
  boundary condition, so their potential stops at the sheath entrance.
- **Self-sputtering.** W ions accelerated through the sheath hit the wall
  with about 3 Z Te + 2 Ti. At high Te or high charge, that is enough for
  W on W to sputter more W (the W on W yield passes 1 near a keV), so
  redeposited W is itself a source. Flan now only counts redeposited atoms.

## What Flan has now

- `Sputtering::Surface` and its elements: the wall surfaces, with a unit
  normal into the plasma and an area, and the plasma cell next to each.
- `Sputtering::Source`: the D+ flux onto each element, the yield, and the
  sputtered atoms' birth (Thompson energy, cosine direction).
- `Sputtering::Deposition`: where atoms come back, how long after being
  sputtered, and their charge. It doesn't know whether a return was prompt.
- `Sputtering::sheath_energy`: the floating sheath plus presheath drop,
  used only for the D+ impact energy. Nothing moves the impurities through
  a sheath.
- The Boris push takes E from the background, interpolated to the particle.
  The time step is fixed (`imp_time_step`), usually 1e-9 to 1e-8 s.
- The slab background (`bkg_source = "slab"`) for regression cases. It
  has B along +X and no wall, and its z is a length, so `Sputtering::Surface`
  doesn't work on it (it needs the toroidal angle).
- `python/rustbca_yields.py` already makes W on W tables
  (`python rustbca_yields.py W W`), but `Sputtering::Source` only accepts
  tables whose projectile is singly charged (D+).

## The physics to add

### Sheath potential

The model is one-dimensional along the wall's normal, with d the distance
from the wall. Two layers sit between the plasma and the wall at an
oblique field:

- **Magnetic presheath** (Chodura layer). Its thickness is a few ion
  sound gyroradii, ρ_s = c_s / Ω_i. Ions are turned from along B toward
  the wall, so the potential drop across it grows as the field gets more
  grazing. With Boltzmann electrons and the ions' normal velocity going
  from c_s sin α to c_s (α the angle between B and the wall), the drop is
  about (Te / e) |ln sin α| (Stangeby, ch. 2).
- **Debye sheath.** Its thickness is a few Debye lengths, λ_D. It takes
  the rest of the floating drop.

The two drops add up to the floating sheath drop in `sheath_energy`, and
the potential is a sum of two exponentials, as in Brooks (1990):

    φ(d) = Δφ_DS exp(-d / L_DS) + Δφ_MPS exp(-d / L_MPS)
    E(d) = (Δφ_DS / L_DS) exp(-d / L_DS) n + (Δφ_MPS / L_MPS) exp(-d / L_MPS) n

with n the unit normal into the plasma. The signs make E point toward
the wall, which pulls positive ions in. L_DS ∝ λ_D and L_MPS ∝ ρ_s.

The shape is written once, and how the drop is split between the two
layers (and the two scale lengths) is a separate, swappable piece:

- **`"brooks"`** (default). Brooks' form, with the Debye fraction
  f_D(α) of the drop from his fit to kinetic results:
  φ(d) = f_D φ₀ exp(−d / 2λ_D) + (1 − f_D) φ₀ exp(−d / ρ_i), with φ₀ the
  floating drop. f_D is near 1 for a field near normal to the wall and
  falls toward 0 at grazing angles. GITR and ERO use this form, so
  Flan's results compare directly with theirs.
- **`"stangeby"`**, for testing how sensitive results are to the split.
  The same shape, with the magnetic presheath drop from fluid theory,
  (Te / e) |ln sin α|, capped at the whole drop, and the rest in the Debye
  sheath. It depends on Ti / Te where the Brooks fit doesn't.

The interface (drops and scale lengths per element and frame) leaves room
for tabulated kinetic (PIC) profiles later, without changing the rest.
The exact fit, the scale lengths, and what happens at very grazing angles
(≲ 1–2°, where the two-layer picture itself is debated, and where the
fluid estimate exceeds the whole drop) are taken from the literature in
Phase S0 and written down with their sources before any code. The field is
added to the background's E only near a wall. The fluid potential
already holds the collisional presheath, so only the two layers above are
added.

The scales are far below a cell (λ_D ~ 10 μm, ρ_s ~ 0.5 mm, cells ~ mm
to cm). So the sheath is an analytic field inside the cell next to the
wall, not something put on the grid.

### Time step near the wall

A W ion crosses a Debye sheath in about 1 ns, which is the usual time
step. Particles close to a wall (d less than a few L_MPS) are sub-stepped,
so that one sub-step moves them a small fraction of the local field's
scale length (L_DS near the wall, L_MPS farther out) and the energy gained
in the sheath is right whatever the global step. This is the same idea as
the collision sub-cycling Flan already does, applied to the push.

### Impact energy and angle

With the sheath on, an ion's energy and angle when it hits the wall come
from its own velocity, since the push has accelerated it. With the sheath
off, the impact energy for self-sputtering is the ion's kinetic energy
plus Z times the sheath drop, the same assumption as for D+.

### Prompt redeposition: Fussmann's model

Fussmann et al. (1995) relate prompt redeposition to
p = λ_iz / ρ, the ionization length of the sputtered atom over the
gyroradius of the ion it becomes. Since both go as the atom's speed,
p = Ω τ_iz, the gyrofrequency of the singly charged ion times the
atom's ionization time, 1 / (n_e ⟨σv⟩_iz), with no dependence on the
emission speed. The fraction redeposited within the first gyration is
about 1 / (1 + p²). The model has no electric field, so it is the limit
Flan must reproduce with the sheath off. The geometry it assumes (field
angle, emission direction and energy) and the exact definition of a
prompt return are fixed in Phase S0 from the paper, and the checker
evaluates the analytic result for the same geometry as the run.

### Self-sputtering

When an impurity ion hits a wall element, it is redeposited (as now) and
sputters Y(E, θ) atoms, from a W on W table, at its impact energy and
angle. Rather than adding particles mid-run, the particle carries on as
one of the atoms it sputtered: it is reborn as a neutral at the same
place on the element, with a Thompson energy (Λ = 1 for W on W) and a
cosine direction, and its weight becomes w Y. Its weight stays the atoms/s
it stands for, so the tallies stay rates:

- Redeposited on the element: w (as now).
- Self-sputtered off the element: w Y.
- Net erosion: D+ source + self-sputtered − redeposited.

When w Y is small, Russian roulette keeps the cost down: the particle
survives with probability w Y / w_min with weight w_min, which keeps the
expected weight. When Y > 1 the weight grows. If the self-sputtering
multiplication (the share coming back times the mean yield) reaches 1,
the chain doesn't end (the self-sputtering avalanche); the run stops with
a message saying so instead of running away.

Reflection of W off W (also from RustBCA) would carry on the same way, as
an atom with the reflected energy. It's left for after the first version
(decision 5).

## Phases

**Phase S0: literature and conventions.** No code. Write down, with
sources, in a draft section of [Sputtering](../physics/sputtering.md):

- Fussmann's model: its assumptions (field angle, emission), the
  definition of p and of a prompt return, and the formula, re-derived.
  `regression/fussmann.py` evaluates it for given B, ne, Te and emission
  distributions, using the same OpenADAS ionization rate Flan uses.
- The sheath model: Brooks' f_D(α) fit and its range of validity, the
  scale lengths (2λ_D and which gyroradius), the Stangeby split as a
  function of α and Ti / Te, the limits (α → 90°, α → 0), and how both
  reduce to the floating drop in `sheath_energy`. Brooks (1990), Chodura
  (1982), Stangeby (2000), and how ERO and GITR set it up, including what
  they do at angles below the fit's range.
- Published prompt-redeposition results with a sheath (e.g., ERO and GITR
  studies, and fits to them) that can serve as reference curves, since
  there is no closed form with a sheath.

**Phase S1: a slab with a wall.** So everything after this can be checked
on the slab, as CLAUDE.md asks of new physics.

- Slab options for a material wall: which face is a wall (`face_wall`),
  and the angle between B and that wall (B tilted in the plane of the
  normal and one in-plane direction). All from `slab_*` options, nothing
  hard-coded.
- `Sputtering::Surface` for straight (Cartesian) zones: the same
  elements with a Cartesian normal and area instead of (R, Z, φ).
- Tests: the surface of the slab wall (normals, areas), the field angle.

**Phase S2: prompt redeposition without a sheath (Fussmann).** Checks
the current code before changing it.

- A tally of prompt returns. Each slot records the time it first
  ionized. A return is prompt if it comes within one gyroperiod of the
  first ionization (or per the definition fixed in S0), and
  `Deposition` counts prompt returns separately, by charge.
- Test-only birth options so the run can match the model's assumptions
  if those differ from the cosine/Thompson birth (decision 6).
- `regression/reg_prompt_redep.cpp` with variants sweeping p from about
  0.1 to 10 through `slab_ne` (and `slab_bmag`), collisions, flow,
  E and the thermal force off. `regression/check_redeposition.py`
  compares the prompt fraction with Fussmann's at each p, within the
  statistical error plus a tolerance set in S0.

**Phase S3: the sheath field (CPU).**

- `src/sheath.cpp`, `include/sheath.h`: the profile φ(d), E(d), and the
  per-element, per-frame parameters (drops, L_DS, L_MPS, α from b · n),
  computed once like `Sputtering::Source`.
- Finding the wall near a particle: each cell next to a wall element
  lists its elements, and d = (x − x_face) · n for the nearest. HERMES
  targets are planes, so d is exact. SOLEDGE3X uses the staircase face the
  particle would be absorbed at, with the element's normal (the
  normal/face mismatch is the same artifact the source already has).
- The extra E added in `update_velocity_cpu` for particles within a few
  L_MPS of a wall.
- Option `sheath_model` (`"off"`, `"brooks"`, `"stangeby"`; `"off"`
  by default until S5 validates it), in all the files listed in "Adding an
  input option".
- Unit tests: ∫E · n dd equals the total drop for both splits; each
  split against its S0 formula over α and Ti / Te; the field points to the wall; zero far
  from it; the lookup finds the right element across zone faces.

**Phase S4: sub-stepping near the wall.**

- Near-wall particles take Boris sub-steps; the absorption itself is still
  done by the boundary check.
- Test: an ion of charge Z started at rest a few L_MPS from the wall, with
  no collisions or ionization, hits it with Z Δφ_total to 1e-3, at time
  steps from 1e-11 to 1e-8 s (energy conservation in a static field).
- Check the cost on a HERMES case: only particles near the wall should
  pay for it.

**Phase S5: prompt redeposition with the sheath.**

- `reg_prompt_redep` variants with the sheath on, same p sweep and a scan
  of field angle and Te, for both `"brooks"` and `"stangeby"`. How much
  the prompt fraction changes between them is the sensitivity to the
  split.
- Once these pass, `sheath_model` defaults to `"brooks"` (the physics is
  then on by default, as for the inductive E∥).
- `check_redeposition.py`: exact checks (ion impact energy against
  Z Δφ + kinetic; Fussmann recovered as the sheath drop goes to 0), and
  comparison with the S0 reference curves (with a looser, stated
  tolerance, since those models differ in detail). The angle of impact
  against the field angle, for the trend the MPS should give (ions turned
  toward the normal).
- Then a HERMES run with and without the sheath, for the change in the
  redeposited share and the W leaving the divertor.

**Phase S6: self-sputtering.**

- A W on W table in `data/sputtering/` from `rustbca_yields.py`, with the
  energy range extended to what 3 Z Te + 2 Ti reaches.
- A second `YieldTable` for the impurity on itself (the projectile check
  in `Source` stays for the D+ table only), impact energy and angle at
  absorption, rebirth with weight w Y, Russian roulette, the avalanche
  stop, a generation counter per slot.
- Options: `sput_self` (`"off"` default), `sput_self_table_file`,
  `sput_self_min_weight`.
- Output: self-sputtered rate per element, net erosion including it,
  weight per generation.
- Unit tests: roulette keeps the expected weight; rebirth puts the atom on
  the element moving into the plasma; the avalanche stop.
- Regression (slab, constant self-yield Y₀ from a synthetic table): the
  weight in generation g+1 over generation g is f Y₀, with f the measured
  return fraction, and the total self-sputtered rate is
  S f Y₀ / (1 − f Y₀) for a D+ source S. Then the real table at a few
  Te, against Y averaged over the impact energies 3 Z Te + 2 Ti and the
  charges the ions land with.

**Phase S7: docs and a full case.** [Sputtering](../physics/sputtering.md)
gets the sheath and self-sputtering sections in place of the "work in
progress" note and the "sheath is still missing" caveat;
[Input Options](../input-options.md) the new options; the regression
README the new cases. A HERMES case with sheath and self-sputtering
shows the net erosion along the targets.

**Later**: the GPU path (sputtering is CPU only today, so this goes with
moving the sputtered source to the GPU), reflection, and a sheath for D+
impact energies from the same model.

## Decisions to make

All decided as recommended.

1. **Sheath model.** Decided: Brooks' two-layer model (Debye sheath plus
   magnetic presheath, with the drop split by his fit to the field angle)
   as the default, as in ERO and GITR, and the Stangeby/Chodura fluid
   split as a second option for sensitivity. Debye sheath only was
   rejected: at grazing angles, where divertor targets are, it puts the
   drop within tens of microns of the wall, where few W atoms ionize, so
   it would barely change prompt redeposition. Tabulated PIC profiles can
   be added later through the same interface.
2. **Geometry order.** Recommended: slab first (S1 to S5), then HERMES
   targets (planes, exact d), then SOLEDGE3X (staircase).
3. **What counts as prompt.** Recommended: a return within one gyroperiod
   of the first ionization, matching Fussmann, settled in S0. Landing
   times and charges stay tallied as now.
4. **Self-sputtering mechanism.** Recommended: weight continuation (the
   particle is reborn as a sputtered atom) over running generations as
   separate passes, since it fits the slots and keeps time order.
5. **Reflection.** Recommended: later, once self-sputtering works.
6. **Test-only birth options** (fixed energy, fixed direction) for the
   Fussmann regression. Recommended only if S0 shows the model's
   assumptions can't be met with the normal cosine/Thompson birth, by
   folding those distributions into the analytic result instead.

## References

- G. Fussmann et al., "High-Z elements as target materials in fusion
  devices", in *Plasma Physics and Controlled Nuclear Fusion Research
  1994*, IAEA, Vol. 2, p. 143 (1995).
- J. N. Brooks, "Near-surface sputtered particle transport for an oblique
  incidence magnetic field plasma", *Phys. Fluids B* **2**, 1858 (1990).
- R. Chodura, "Plasma–wall transition in an oblique magnetic field",
  *Phys. Fluids* **25**, 1628 (1982).
- P. C. Stangeby, *The Plasma Boundary of Magnetic Fusion Devices*, IOP
  Publishing (2000).
