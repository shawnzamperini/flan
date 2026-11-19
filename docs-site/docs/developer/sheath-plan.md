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
- Option `sheath_model` (`"off"`, `"brooks"`, `"borodkina"`; `"off"`
  by default until S5 validates it), in all the files listed in "Adding an
  input option".
- Unit tests: ∫E · n dd equals the total drop for both models; f_D
  against its S0 formula over θ and Ti / Te, with the floor past the
  critical angle; `"borodkina"` continuous in φ and E at ξ_mps and against
  values from its equations; the field points to the wall and is zero far
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
