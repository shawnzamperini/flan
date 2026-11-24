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
- The slab background (`bkg_source = "slab"`) for regression cases. Before
  S1 it had B along +X and no wall, and its z is a length, so
  `Sputtering::Surface` didn't work on it (it needed the toroidal angle).
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

Both models split the floating drop the same way, with the Debye
fraction f_D = 1 − ln cos θ / λ_w (θ the field angle from the normal,
λ_w = eφ_fl / Te; Stangeby 2012, Borodkina and Tsvetkov 2015). Past
Stangeby's critical angle this would be negative, so f_D is held at a
small floor there, set from Coulette and Manfredi's kinetic results,
which show a smooth change rather than a sharp cut-off. The models differ
in the shape of φ(d), behind one interface that gives φ and E at a
distance d from an element:

- **`"brooks"`** (default). Brooks' form,
  φ(d) = f_D φ₀ exp(−d / 2λ_D) + (1 − f_D) φ₀ exp(−d / ρ_i), with φ₀ the
  floating drop. GITR and ERO use this form, so Flan's results compare
  directly with theirs. Brooks used a constant f_D and GITR an f_D curve
  with no published source; Flan uses the split above
  ([S0 findings](#phase-s0-findings)). The gyroradius definition (Brooks'
  √(kT / m) / Ω or GITR's √(2kT / m) / Ω) is fixed in S3.
- **`"borodkina"`**, for testing how sensitive results are to the shape
  of the profile. Borodkina and Tsvetkov's analytic profile: a Debye layer
  solved from Poisson's equation next to the wall, then an exponential
  magnetic presheath with decay length ρ_s sin θ. Past the critical angle
  it is presheath only, anchored at the wall. It has been checked against
  PIC (SPICE2), and GITR benchmarks against it.

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

**Phase S0: literature and conventions.** Mostly done; see
[Phase S0 findings](#phase-s0-findings) for what was found and what's
left. No code apart from `regression/fussmann.py`. Write down, with
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
on the slab, as CLAUDE.md asks of new physics. Done: `slab_wall` makes one
face of the straight slab a wall (`face_wall`), and `slab_b_angle` turns B
from +X toward −Y, so with the wall at y min it's the angle between B and
the wall. Zones can be marked Cartesian (`Zones::set_cartesian`), and on
those `Surface::find_wall_faces` takes the corners from the zone's grid,
with R, Z and phi holding X, Y and Z (`SurfaceElement::cartesian`,
`normal_xyz`). Tests in `test_sputtering_cpu.cpp` (`SlabWall`): the
surface at y min and x max, the field angle and flow, the flux
n u sin α onto the wall and the births, deposition on the element under
the atom, and the option checks. Wall faces on the slab are CPU only.

- Slab options for a material wall: which face is a wall (`face_wall`),
  and the angle between B and that wall (B tilted in the plane of the
  normal and one in-plane direction). All from `slab_*` options, nothing
  hard-coded.
- `Sputtering::Surface` for straight (Cartesian) zones: the same
  elements with a Cartesian normal and area instead of (R, Z, φ).
- Tests: the surface of the slab wall (normals, areas), the field angle.

**Phase S2: prompt redeposition without a sheath (Fussmann).** Checks
the current code before changing it. Done:

- `Slots::prompt_by`: when an atom first ionizes,
  `OpenADAS::ioniz_recomb` sets it to that time plus 2π m / (eB) at the
  local field. `Deposition` counts a return before then as prompt
  (`prompt`, `prompt_charge`, `prompt_total` in the output).
- The model's assumptions can't all be met with the normal physics: W+
  ionizes again about 2π (S₁/S₀)/p times per gyration. So decision 6's
  test-only option was needed after all, but for the charge rather than
  the birth: `imp_iz_max_charge` (CPU only) stops ionization at W+. The
  cosine/Thompson birth is fine, since p doesn't depend on the speed.
- B is at 0.01° to the wall, so the plasma reaches it (the source needs
  a flux) while an ion moves only ~1e-3 ρ along B per gyration.
- `reg_prompt_redep` (p01 to p10, b4, and `full` without the charge
  limit) and `check_redeposition.py`, which takes p from the saved
  inputs and the same OpenADAS table, interpolated as Flan does.
  Tolerance: 4 standard errors plus 0.01.
- Unit tests: `IonizRecomb.PromptBy`, `IonizRecomb.MaxCharge`,
  `SlabWall.PromptTally`.

The plan below was:

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

**Phase S3: the sheath field (CPU).** Done:

- `Sheath::Sheath` (`sheath.h`, `sheath.cpp`) holds `Params` for every
  wall element and time frame (Te, Ti, n_e and B of the element's plasma
  cell; θ from B and the element's normal), and for each plasma cell next
  to the wall, its elements. `potential`, `field` and `range` give φ(d),
  E · n(d) and how far the field reaches for either model.
- The gyroradius is ρ_s = c_s / Ω_i with c_s = √(e(Te + Ti)/m_i): Borodkina's
  ρ_iCs, and GITR's √(2kT/m)/Ω at Ti = Te (Brooks' R_DT is √2 smaller).
- No floor on f_D past the critical angle: f_D goes to 0 continuously
  there (only its slope jumps), and the kinetic results give no value to
  set a floor from. Past it both models put the whole drop in the
  presheath. Easy to add as an option if S5 shows it matters.
- `"borodkina"` at normal incidence: Q = −λ_w and ξ_mps is infinite (one
  Debye layer, decaying over 1/a ≈ 2 λ_D, as in Brooks). φ and E are
  continuous at ξ_mps at every angle checked (10° to 85°). Within about a
  degree of the critical angle a and Q take odd values, but the Debye
  layer holds only a few percent of the drop there.
- The distance from an element is measured across its plasma cell: the
  fraction of the way from the element's face to the opposite face, times
  the cell's thickness along the element's normal (2 (c − m) · n). Exact
  on the slab and on rectangular cells; on the staircase it is the
  distance along the outline's normal. The field is only added in the
  cells next to the wall, so it is cut off at the far side of a cell
  thinner than the sheath; `Sheath` counts such elements and
  `print_summary` warns.
- `update_velocity_cpu` adds `Sheath::field_at` to the interpolated E.
  `flan.cpp` builds the `Sheath` from the sputtered source's surface (or
  its own) when `sheath_model` isn't `"off"`.
- `sheath_model` defaults to `"brooks"` (the user's call, ahead of S5).
  Backgrounds without walls (`Sputtering::has_surfaces`) go without it,
  and a GPU run goes without it with a warning. `reg_prompt_redep` turns
  it off, since Fussmann's model has no sheath.
- Tests: `test_sheath_cpu.cpp` (the scales, f_D over θ and the critical
  angle, ∫E · n = φ_wall and E = −dφ/dd for both models over 0° to 89°,
  Brooks' and Borodkina's forms, the field on the slab wall, and an ion
  falling through the field gaining φ(start) − φ_wall), and
  `SheathField.DiagonalWall` (toroidal 45° staircase: distance, θ, and
  the normal turning with φ).

The plan was:

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
- Option `sheath_model` (`"off"`, `"brooks"`, `"borodkina"`; `"off"`
  by default until S5 validates it), in all the files listed in "Adding an
  input option".
- Unit tests: ∫E · n dd equals the total drop for both models; f_D
  against its S0 formula over θ and Ti / Te, with the floor past the
  critical angle; `"borodkina"` continuous in φ and E at ξ_mps and against
  values from its equations; the field points to the wall and is zero far
  from it; the lookup finds the right element across zone faces.

**Phase S4: sub-stepping near the wall.** Done:

- `Sheath::substeps` gives each charged particle that can reach the
  field this step (d − v_max dt within 7 scale lengths) a number of
  sub-steps N so none moves it more than 2% (`substep_frac`) of the local
  scale length: the Debye layer's if it can reach that layer this step,
  the presheath's if not. v_max is its speed plus what the rest of the
  drop could add. At least 2, at most 20000 (`max_substeps`).
- `Boris::substep_cpu` pushes those particles in N sub-steps of dt/N
  (Boris update, move, find the cell). The leapfrog's half-step offset is
  moved from dt/2 to h/2 at the start and back at the end; Boris is
  reversible, so these cancel on consecutive near-wall steps. A particle
  leaving the grid gets a kick of (f − ½)h at its last position inside (f
  the fraction of the sub-step at which it crossed), so it hits the wall
  with its velocity at the crossing. `update_velocity_cpu` and
  `step_cpu` skip the sub-stepped particles' push and move; `step_cpu`
  still advances their time.
- The Boris loop body is now `Boris::push_one_cpu`, used by both.
- Test `SheathSubstep.ImpactEnergy`: a W ion of charge 1 or 4 released at
  rest 4 scale lengths out, B at 30° to the wall, hits it with
  Z (φ(start) − φ_wall) to 1e-3 for both models at dt = 1e-10, 1e-9 and
  1e-8 s (1e-11 s passes too, but takes over a minute in a debug build).
  Without the crossing kick the error was up to 1%.
  `SheathSubstep.NeededAtLargeSteps` shows that the plain push is off by
  more than 5% at 1e-8 s.
- Cost: on the slab prompt redeposition case, where nearly every particle
  is near the wall, the run takes about 1.5 times as long. On a HERMES
  turbulence background with a sputtered target source, the push's time
  is unchanged (the run was even shorter, as more atoms come back
  promptly).
- On that background about half the target elements have a plasma cell
  thinner than 7 presheath lengths, leaving out up to about 4% of the
  drop past the cell's far side. `print_summary` now reports the largest
  such fraction. Carrying the field into the next cell out would fix it;
  left for now.

The plan was:


- Near-wall particles take Boris sub-steps; the absorption itself is still
  done by the boundary check.
- Test: an ion of charge Z started at rest a few L_MPS from the wall, with
  no collisions or ionization, hits it with Z Δφ_total to 1e-3, at time
  steps from 1e-11 to 1e-8 s (energy conservation in a static field).
- Check the cost on a HERMES case: only particles near the wall should
  pay for it.

**Phase S5: prompt redeposition with the sheath.** Done:

- Impact tallies in `Deposition` (`Sputtering::Impact`, passed from the
  boundary check): the energy and the angle from the element's normal of
  each absorbed particle, as histograms, means and the energy on each
  element (for S6), and the mean time since birth. Saved in the
  sputtering group and printed in the summary.
- Options `sheath_fd` and `sheath_drop` fix f_D and the drop (over Te),
  for Brooks' reference case and for checks with a vanishing drop.
  `"borodkina"` uses −Δn λ_mps = (λ_mps / L)², so a given f_D carries
  through to it too.
- `reg_sheath_angles` (a45, a60) and `check_sheath_angles.py`: D⁺ test
  ions from 10 ρ_s through the `"borodkina"` field, pushed as in the main
  loop. Energy is conserved to 4e-4 of the drop or better for every ion.
  The impact angle distribution matches Borodkina and Tsvetkov's at 45°
  (mean 30.4° against 32.4°) and has their shape and width at 60° but
  about 4° nearer the normal (31.5° against 35.9°). The time step and the
  profile (`"brooks"` gives 31.1°) don't change that; the start does (the
  paper only says "a shifted Maxwell velocity distribution"; shifting by
  0.5 or 2 c_s instead of c_s gives 29° or 39°), so the mean's tolerance
  is 5°.
- `reg_prompt_redep` variants: `drop0` (sheath on, drop 1e-4 Te) recovers
  Fussmann and is checked; `<model>_a<angle>` at 2° to 20° are printed and
  plotted. At p ≈ 1 the sheath raises the prompt fraction from 0.59 to
  0.72 to 0.86, less as B gets steeper, and the two profiles agree to
  about 0.015.
- `brooks_ref` is printed beside Brooks' Table 1, not checked: Flan gives
  charge 1.4, 105 eV, 13° from the normal and 0.04 μs, against his 2.6,
  553 eV, 72° and 0.6 μs. Flan's W⁰ ionization length there is about
  0.06 mm (OpenADAS) against his 0.40 mm mean ionization height, so its W
  ionizes deep in the presheath. The rates differ, not the sheath.
- HERMES (TCV-X21 turbulence background, sputtered targets, 1000 W for
  20 μs): prompt redeposition 53% without the sheath and 67% with
  `"brooks"`. The W leaving the divertor wasn't measured.
- `sheath_model` was already made `"brooks"` by default before S5.

The plan was:


- `reg_prompt_redep` variants with the sheath on, same p sweep and a scan
  of field angle and Te, for both `"brooks"` and `"borodkina"`. How much
  the prompt fraction changes between them is the sensitivity to the
  shape of the profile.
- D⁺ test ions pushed through the `"borodkina"` field reproduce
  Borodkina and Tsvetkov's impact-angle distributions (Fig. 2a; 45° and
  60°), as a check on the field and the near-wall push.
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

**Phase S6: self-sputtering.** Done:

- `Sputtering::SelfSource` (CPU only). `Deposition::add_at_wall_face`
  and `add_into_wall_cell` return the element a particle landed on, and
  the boundary check passes it to `SelfSource::rebirth` (reached through
  `Deposition::self`, so nothing else is threaded through the main
  loop).
- Impact energy and angle (`self_impact`): from the velocity with the
  sheath field on; without it, an ion gets q φ_wall added to its motion
  along the normal, with φ_wall the floating drop alone (as in the sheath
  model, or `sheath_drop` Te), not the D⁺ 0.5 Te presheath. Neutrals keep
  their own energy. `Source::phi_wall` holds the drop per frame and
  element.
- Rebirth: Y from the table, the element's self-sputtered rate gets w Y,
  and the weight w Y goes through roulette (`roulette`, below
  `sput_self_min_weight` times the source particles' weight). A survivor
  becomes a neutral where it hit, a thousandth of a cell into the plasma
  cell, with a Thompson energy cut off at Λ E (Λ = 1 for W on W) and a
  cosine direction. Its time since birth and `prompt_by` are reset; the
  background time goes on. `Slots::generation` counts its generations.
- Avalanche: a particle that would pass `sput_self_max_gen` (1000) stops
  the run after the boundary check, with the multiplication measured
  along the chains. With roulette a chain survives each generation with
  probability about min(f, f Y), so long chains only happen when that is
  near 1.
- Output (sputtering group): `self_sputtered` and `net_erosion_self` per
  element, `generation_weight` and `generation_landed`, the
  multiplication T / (S + T), the chain estimate, mean yield, impact
  energy and angle, and what roulette ended. `fates_self_sputtered`
  counts the rebirths.
- Tests: `SelfSputtering.Roulette`, `SelfSputtering.Impact`,
  `SputteringTable.SelfRebirth` (place, direction, weight, tallies,
  sheath on and off, below threshold, generation cap),
  `SelfRoulette`, `SelfChecks`, `SelfBoundary` (through
  `check_bounds_cpu`).
- `reg_self_sputter` y05 and y125 (constant Y0 from
  `make_self_tables.py`) and `check_self_sputter.py`: W_{g+1} = Y0 L_g
  (exact without roulette, binomial with it), the tallies agree,
  f₁ = f₀ (0.585 against 0.596 at Y0 = 0.5, 0.583 against 0.588 at 1.25),
  and T = S f₀ Y0 / (1 − f₁ Y0) (T / S = 0.416 and 2.749). All pass.
  About 4 minutes each on 16 cores.
- `data/sputtering/W_on_W.h5` (RustBCA 3.0.0, 10 eV to 10 keV) and the
  `w_te<Te>` variants (full ionization, "brooks" sheath, B at 2 degrees),
  printed by the checker. At Te = 20 eV the mean yield of what lands is
  0.012 and T / S = 0.009 (3 generations); at 50 eV it is 0.091 and
  0.089 (7 generations). The checker's rough estimate, every landing
  particle an ion at q φ_fl + 2 Ti, gives 0.032 and 0.31: the yield is
  steep there, and the mean impact energy (58 and 124 eV) is well below
  that estimate's. About 3 and 2 minutes on 16 cores.

The plan was:

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

**Phase S7: docs and a full case.** Done:

- [Sputtering](../physics/sputtering.md): "The sheath" and
  "Self-sputtering" sections, the impact tallies under Redeposition, the
  W on W table, the sheath and self-sputtering regressions under
  Validation, and the "work in progress" note and the "sheath is still
  missing" caveats replaced. [Input Options](../input-options.md) and the
  regression README got the new options and cases in S3 to S6.
  [Code Structure](code-structure.md) says where the source, the sheath
  and self-sputtering hook into the main loop.
- `python/plot_sputtering.py` prints the impact means, the self-sputtered
  rate, multiplication and generations, counts self-sputtering in the
  erosion maps, and plots the sputtered, self-sputtered, redeposited and
  net rates along each target (`_targets.png`, for runs with at most 8
  targets).
- HERMES (TCV-X21 turbulence background, sputtered targets, `"brooks"`
  sheath, 10⁴ W for 20 μs, W_on_W.h5), with self-sputtering off and on:
  T / S = 0.018 (mean yield 0.026 at 75 eV and 18° from the normal, down
  to generation 3); 68% of what is sputtered comes back either way; the
  net erosion over the targets goes from 3.16e15 to 3.22e15 atoms/s, about
  2%, about the statistical noise at 10⁴ particles. Self-sputtering
  is strongest where the D sputtering is. 14 and 12 minutes on 16 cores.

The plan was: [Sputtering](../physics/sputtering.md)
gets the sheath and self-sputtering sections in place of the "work in
progress" note and the "sheath is still missing" caveat;
[Input Options](../input-options.md) the new options; the regression
README the new cases. A HERMES case with sheath and self-sputtering
shows the net erosion along the targets.

**Phase S8: D⁺ impacts from the sheath model.** The sputtered source
takes the D⁺ impact energy from a formula (the floating drop, 0.5 Te of
presheath and a thermal 2 Ti) and every ion hits at `sput_angle_val`
(60° by default). Pushed through the sheath field, D⁺ hits with about
the same mean energy, but nearer the normal (about 30° in
`reg_sheath_angles`), with energy and angle correlated, and near
threshold the D on W yield depends on the angle. So the source should
average the yield over the energies and angles the sheath model gives.

The 1-D problem, in units of ρ_s, c_s and Te, depends only on the field
angle θ_B (which also fixes f_D), Ti / Te and λ_D / ρ_s. So test ions are
pushed for a grid of those, not per element and frame, and each element
and frame uses its grid point's samples.

- **D1. Test-ion engine.** Done:
  - `Sheath::sample_impacts` pushes D⁺ test ions (charge 1) in the 1-D
    sheath frame: wall at y = 0, uniform B at θ from the normal pointing
    at the wall, E = `field(p, model, y)` along the normal. They start at
    10 ρ_s from a Maxwellian at Ti shifted by c_s along B, keeping those
    moving toward the wall along B, as in `reg_sheath_angles`. It returns
    each hit's energy, angle from the normal, and weight v · b at the
    start (the ion's share of the flux into the sheath), with how many
    were lost and the largest energy error. Each ion has its own random
    stream, so results don't depend on threads.
  - The push is the main loop's: the Boris update is now
    `Boris::kick` (shared by `push_one_cpu`, same arithmetic in the same
    order) and the sub-step rule a free `Sheath::substeps(Params, ...)`
    (called by the member), with the same half-step offsets and crossing
    kick as `Boris::substep_cpu`. The step is 1/100 of a gyroperiod; at
    1/50 the a60 energy error was 8.6e-4 of the drop, at 1/100 3.8e-4
    (the main loop's), for the same run time. Ions are lost after 2000
    gyroperiods or going back past twice their start.
  - `Params` now carries Ti, |B| and the ion mass, so the sampler needs
    nothing else and can't be given a Ti that disagrees with ρ_s.
  - Tests (`test_sheath_cpu.cpp`, `SheathImpacts.*`): energy conserved to
    1e-3 of the drop and at most 2% lost for both profiles at 0°, 30°,
    60°, 85° and 89°; at normal incidence the speed across B is the
    starting Rayleigh one; ions hit nearer the normal than the field
    (64° to 68° at 85°); reproducible by seed.
  - `reg_sheath_angles` `a45_sampler` and `a60_sampler`, checked by
    `check_sheath_angles.py` as the main-loop runs are and against them:
    mean angles 30.40° and 31.51° against the main loop's 30.42° and
    31.50°, standard deviations within 0.12°. 2 to 4 seconds against 2
    to 5 minutes. Weighted by the flux, the means are 31.0° and 34.3°.
- **D2. Table and the source.** Done:
  - `Sheath::ImpactTable`: points on a grid of 0.5° in θ_B, 10% in
    Ti / Te and factors of 2 in λ_D / ρ_s (θ_B held at 89.5° at most),
    added only where an element and frame falls (`add`), then sampled
    (`build`) with 2000 ions each and a seed from the point's grid
    indices, so its samples don't depend on the run's other points.
    Each point is sampled at its grid values with Te = 10 eV and B = 1 T
    (n_e set by λ_D / ρ_s); energies are kept over Te and weights
    normalized.
  - `Sputtering::Source` takes the D⁺ impact model as a new constructor
    argument (the formula by default; D3 adds the option). With a sheath
    model, each element and frame gets its point from `make_params`
    (with `sheath_fd` and `sheath_drop`), and Y_eff, the mean impact
    energy and angle, and the share above the table are averages over
    the point's samples weighted by their flux share, with energies times
    the element's Te. Births draw the ion from those samples weighted by
    w Y (O(2000) yield lookups a birth), so `Birth` now carries the
    ion's angle too (`a_impact`). Elements with no field keep the
    formula. Γ is unchanged. `Source::impact_angle` and `impacts`, and
    the summary prints the mean angle and the table.
  - The mean D⁺ energy is Ti higher than the formula's: the flux into
    the sheath of a Maxwellian drifting at c_s carries m c_s² / 2 + 2.5 Ti
    = 0.5 Te + 3 Ti on top of the drop, against the formula's
    0.5 Te + 2 Ti. At Te = Ti = 20 eV, 128 eV against 107 eV (and
    0.5 Te + 3 Ti predicts 126.8 eV).
  - Tests: `SheathImpactTable.Points` and `ScalesWithTe` (a real sheath
    at 40 eV and 3 T with a point's three numbers gives its energies
    times 40 and its angles); `SputteringTable.SheathImpacts` (constant
    yield gives Y₀ and the formula's rate exactly; one point on the slab;
    angle between 5° and 70° with B 80° from the normal; energy within 3%
    of the drop + 0.5 Te + 3 Ti) and `SheathImpactsStep` (a step yield
    gives the flux share above it, and births only pick ions above
    it).
- **D3. Option and output.** `sput_impact_opt` = `"formula"` (now's
  behaviour, the default until D4 passes) or `"sheath"`, which uses
  `sheath_model`'s profile (and `sheath_fd`, `sheath_drop`) and stops
  with an error if it is `"off"`. Output: the mean D⁺ impact angle per
  element and frame, D⁺ energy and angle histograms over the run, the
  grid and sample count; the start-of-run summary prints the mean angle.
- **D4. Regression.** A constant yield gives Y_eff = Y₀ exactly; a yield
  linear in E and flat in angle gives back the mean impact energy
  (within a few percent of the drop + 0.5 Te + 3 Ti, see D2); at normal
  incidence the source matches `"formula"` with `sput_angle_val` = 0;
  the 45° and 60° histograms match Borodkina's. D on W at a few Te and
  field angles is printed, `"formula"` beside `"sheath"`.
- **D5. Docs and a full case.** The D⁺ impact section in
  [Sputtering](../physics/sputtering.md) and the option in
  [Input Options](../input-options.md); TCV-X21 with both options (gross
  source, its spread along the targets, net erosion). Then `"sheath"`
  becomes the default.

Decisions (2026-11-23): the c_s-shifted start, fixed, as in
`reg_sheath_angles`; the table built at startup, not a file; `"sheath"`
follows `sheath_model` and refuses `"off"`; `"formula"` stays the
default until D4 passes; the E × B share of Γ gets the same
distribution, stated as a limitation.

**Later**: the GPU path (sputtering is CPU only today, so this goes with
moving the sputtered source to the GPU) and reflection.

## Phase S0 findings

### Fussmann's model

- **Source.** The derivation isn't in an accessible copy of Fussmann et
  al. (1995). Reimold (2011, IPP 17/27, sections 2.1.1.3 and 3.3.1.2)
  gives the model in full, attributes the derivation to R. Dux (private
  communication), and implements it in EIRENE. Later papers cite Fussmann
  et al. (1995) for it.
- **Assumptions.** The field is parallel to the wall, sputtered atoms leave
  with a cosine distribution, and there's no electric field. Motion along B
  and further ionization during the first gyration are left out. An ion is
  promptly redeposited if its first gyro-orbit reaches the wall.
- **Parameter.** p = λ_iz / ρ = q B / (m n_e S), with S the ionization
  rate coefficient of the atom and ρ the gyroradius of the singly charged
  ion at the atom's full speed. p doesn't depend on the emission energy.
- **Result.** For a fixed ionization length x ρ, r(x) = 1 / (1 + x²).
  Averaged over the exponential distribution of ionization lengths (what
  Flan samples), ⟨r⟩(p) = p⁻¹ ∫₀^∞ e^(−x/p) / (1 + x²) dx. At p = 1 that's
  0.62, not the 0.50 of 1 / (1 + p²). The S2 checker uses ⟨r⟩.
- **Checked.** A Monte Carlo of the geometry reproduces 1 / (1 + x²) to
  the statistics (`regression/fussmann.py`, `geometry_check`), and ⟨r⟩
  reproduces Reimold's analytic value (37% at p = 2.55).
- **For S2.** Flan follows ionization beyond W⁺. The model leaves that out,
  so the regression needs n_e and Te where W⁺ rarely ionizes again within
  one gyration. The checker estimates that fraction and flags the variants
  where it isn't small. It also takes S at the run's Te; the repo has no
  Python OpenADAS reader, so S2 adds one or saves S in the output.
- **Self-sputtering.** Reimold (eq. 2.18) gives the net yield summed over
  the chain, Y_eff = (1 − ⟨r⟩) Y_D / (1 − ⟨r⟩ Y_W), which runs away when
  ⟨r⟩ Y_W ≥ 1. It goes into the S6 checks (`fussmann.y_eff`).

### Brooks (1990)

Read from the ANL report version (ANL/FPP/TM-244, 1989), which is the
same paper as Phys. Fluids B 2, 1858 (1990).

- **Potential.** φ(z) = φ₁ exp(−z / 2λ_D) + φ₂ exp(−z / R_DT), with
  φ₁ + φ₂ = φ₀, λ_D the Debye length at the sheath entrance, and R_DT the
  D-T ion gyroradius (eq. 11). E has only a normal component, E_z = −dφ/dz.
  Electrons are Boltzmann, Ti = Te, and the flow before the sheath is
  sonic.
- **Split.** f_D = 1 − φ(6λ_D) / φ₀ is the fraction of the drop in the
  Debye sheath, whose width is taken as 6 λ_D (eq. 12). f_D is a constant
  input: 0.25 in the reference case, "typical of results for highly
  oblique field angles" (from Chodura). Brooks has no f_D(angle) fit.
- **Angles.** The field angle ψ is measured from the surface normal; the
  reference case has ψ = 87°.
- **Gyroradius.** R_DT = √(kT / m_DT) / Ω_DT. This reproduces his
  "3 R_DT = 0.84 mm at Te = 75 eV, B = 5 T".
- **Reference case for S5.** B = 5 T, ψ = 87°, eφ₀ = −3 kTe, f_D = 0.25,
  Te = Ti = 75 eV, n_e = 10²⁰ m⁻³, D-T, collisions on (Table 1). For
  redeposited W he gets:

  | Quantity | Value |
  |---|---|
  | Mean ionization height | 0.40 mm |
  | Transit time | 0.60 μs |
  | Elevation angle | 18° |
  | Charge | 2.6 |
  | Impact energy | 553 eV |

  Figs. 4 to 8 give the mean charge, energy and self-sputtering yield for
  Te from 25 to 150 eV and n_e from 0.25 to 4 × 10²⁰ m⁻³, and how they
  change with f_D and φ₀. His ionization rates and charge limit (K ≤ 11)
  differ from Flan's OpenADAS, so expect trends, not exact numbers.

### GITR's form and f_D(α)

- **Form.** GITR (Kumar et al., arXiv:2609.09408, sec. 4.3 and app. B)
  uses E_n = f_D φ₀ / (2λ_D) e^(−d/2λ_D) + (1 − f_D) φ₀ / ρ_i e^(−d/ρ_i),
  with θ_B measured from the normal.
- **Gyroradius.** Its ρ is √(2 kTi / m_i) / Ω_i, √2 larger than Brooks'
  R_DT. This reproduces their λ_D = 1.919e-5 m and ρ = 2.229e-4 m at
  Te = Ti = 20 eV, n_e = 3 × 10¹⁸ m⁻³, B = 4.1 T.
- **f_D(θ_B).** Their Fig. B1 inset reads about 1 at 0°, 0.8 at 60°,
  0.45 at 80° and 0.1 at 89°. The paper gives no formula or source for
  it, beyond citing Brooks (1990), who has none. They benchmark the field
  against Borodkina et al. (2015, 2016): close for θ_B ≤ 60°, small
  differences at 80° and above.

### Stangeby's critical angle and the kinetic results

- **Fluid result.** Stangeby (Nucl. Fusion 52, 083012, 2012) finds that
  below a critical angle there is no Debye sheath: the flow normal to the
  wall stays subsonic, and the whole drop is in the Chodura layer. With α
  measured from the wall, sin α_c = c_s / ⟨v∥e⟩_W, which is
  √(2π (m_e / m_i)(1 + Ti / Te)) for half-Maxwellian electrons. That's
  3.4° for D at Ti = Te, and around 5° for ITER-like conditions.
- **Kinetic simulations.** Coulette and Manfredi (arXiv:1509.04479,
  2016) confirm this with kinetic simulations: the space charge at the wall
  fades below about 5°, and the field spreads over several ion Larmor
  radii. They also find no sharp critical angle; the change is smooth, and
  the total drop barely depends on α.
- **The Stangeby split.** The magnetic presheath drop is
  (Te / e) |ln sin α|, and f_D = 1 − |ln sin α| / (e φ_fl / Te). This
  reaches 0 exactly at α_c, since the presheath drop then equals the
  floating drop, so the split has the critical angle built in.
- **Against GITR's curve.** At θ_B = 30°, 60°, 80° and 85° the Stangeby
  split gives 0.95, 0.76, 0.38 and 0.14, close to GITR's curve. Below α_c
  it gives 0, where GITR's stays at about 0.1, nearer the smooth kinetic
  behaviour.

### Borodkina and Tsvetkov (2015)

Physics Procedia 71, 25 (2015), the conference version of the analytic
model that GITR benchmarks against. α is measured from the surface
normal, and distances ξ are in Debye lengths r_d. ρ_s is the Larmor
radius at the ion acoustic speed, c_s / Ω_i with c_s = √(e(Te + Ti) / m_i).
Potentials λ are e(φ − φ₀)/kTe.

- **Wall potential.** λ_w = ½ ln(2π (m_e / M_i)(Te + Ti) / Te). This is
  the floating drop, with no presheath term.
- **Split.** λ_mps = ln cos α, Stangeby's magnetic presheath drop. So
  f_D = 1 − λ_mps / λ_w, the same split as Stangeby's.
- **Magnetic presheath** (eq. 5). The potential is an exponential,
  λ = λ_mps exp(−(ξ − ξ_mps) / L), with
  L = √(−λ_mps / Δn_mps) = ρ_s sin α / r_d. That is, the decay length is
  ρ_s sin α, which grows with the angle. Brooks and GITR use a constant
  gyroradius instead.
- **Debye sheath** (eqs. 6, 8). λ = λ_w + Q (1 − e^(−aξ)), with a and Q
  closed-form functions of λ_w, λ_mps and α.
- **Boundary** (eq. 9). ξ_mps comes from continuity of the potential.
- **Checked.** I evaluated the model for D at Te = Ti = 30 eV,
  n = 10²⁰ m⁻³ and B = 3.2 T (their Fig. 1 case), at 30° to 85°. The
  field is continuous at ξ_mps, which also checks the transcription of
  the scanned equations. The Debye layer is 8 to 13 λ_D wide.
- **Where it stops.** At α = 0 there is no presheath and the formulas
  divide by zero; use the Debye part alone. Past Stangeby's critical
  angle (86.6° here), λ_mps < λ_w and a, Q and ξ_mps go negative. The
  natural extension is no Debye sheath, with the presheath exponential
  anchored at the wall and capped at λ_w.
- **Benchmarks.** The profile agrees with Chodura's and Stangeby's
  solutions and with SPICE2 PIC at 80° (Fig. 1). The impact-angle
  distribution of D⁺ agrees with SPICE2 at 60° and with Kawamura et al.
  PIC at 45° (Fig. 2a). Those are D⁺ test ions, not W, but Flan can
  reproduce them in S5 by pushing D⁺ test ions through the sheath field.
- **GITR's curve.** GITR's f_D stays near 0.1 at 89°, where this model has
  no Debye sheath. So GITR isn't simply 1 − λ_mps / λ_w, and its low-angle
  values come from somewhere unpublished.

### Still to do in S0

- **Where to get the `"brooks"` f_D(α) from.** Decision 7 below.
- **Borodkina et al. (2016)**, Contrib. Plasma Phys. 56, 640. It may give
  more benchmarks, but the model is in the 2015 paper.
- **Reference curves with a sheath.** Published prompt-redeposition
  results (ERO, GITR, Chankin et al. 2014) for S5. Brooks' reference case
  above is the first one.

## Decisions to make

All decided. 1 to 7 as recommended.

1. **Sheath model.** Decided: Brooks' two-layer model (Debye sheath plus
   magnetic presheath) as the default, as in ERO and GITR, and a second
   model for sensitivity (settled in decision 7). Debye sheath only was
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

7. **Where `"brooks"` gets f_D(α), and what the sensitivity option is.**
   Decided: (a). Brooks used a constant, and GITR's curve has no published
   formula.
   Stangeby and Borodkina split the drop the same way,
   f_D = 1 − ln cos α / λ_w (α from the normal). Options:
   (a) `"brooks"`: Brooks' two-exponential shape with that f_D, held at a
   small floor past the critical angle (the kinetic results show no sharp
   cut-off). Replace the planned `"stangeby"` option with `"borodkina"`:
   the same split, but with the Poisson-shaped Debye layer and the
   ρ_s sin α presheath scale. The two then differ only in the shape of the
   profile, which is the uncertainty worth testing;
   (b) as (a), but `"borodkina"` is the default and `"brooks"` the
   comparison with GITR and ERO;
   (c) digitize GITR's curve (Kumar et al., Fig. B1) for f_D.
   Recommended: (a). It is Brooks' form, which GITR and ERO use, with an
   f_D that has published physics behind it. `"borodkina"` is a published,
   PIC-benchmarked alternative with its own S5 test (the D⁺ impact
   angles). (c) copies a curve whose source we can't name.


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
- P. C. Stangeby, "The Chodura sheath for angles of a few degrees between
  the magnetic field and the surface of divertor targets and limiters",
  *Nucl. Fusion* **52**, 083012 (2012).
- D. Coulette and G. Manfredi, "Kinetic simulations of the Chodura and
  Debye sheaths for magnetic fields with grazing incidence", *Plasma Phys.
  Control. Fusion* (2016), arXiv:1509.04479.
- F. Reimold, "Numerical simulations of tungsten impurity transport in
  ASDEX-Upgrade with EIRENE", IPP Report 17/27 (2011).
- J. N. Brooks, ANL/FPP/TM-244 (1989), the report version of Brooks (1990).
- A. Kumar et al., "RF-specific tungsten erosion and global transport in
  ITER under neon seeding", arXiv:2609.09408 (2026).
- I. E. Borodkina and I. V. Tsvetkov, "Modelling of charged particle
  dynamics in the sheath and plasma-facing surface sputtering", *Physics
  Procedia* **71**, 25 (2015).
- I. E. Borodkina, M. Komm, I. V. Tsvetkov, *Russ. Phys. J.* **58**, 438
  (2015); I. Borodkina et al., *Contrib. Plasma Phys.* **56**, 640 (2016).
