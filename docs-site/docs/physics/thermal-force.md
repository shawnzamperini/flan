# Thermal Force

An impurity ion in a plasma with a temperature gradient along the magnetic
field is pushed toward the hotter side. This is the thermal force. For a
heavy impurity of charge state $Z$ in a collisional plasma it is

$$
F_\text{th} = 0.71\,Z^2\,\nabla_\parallel T_e + 2.65\,Z^2\,\nabla_\parallel T_i ,
$$

with the temperatures in energy units. It grows as $Z^2$, like the friction
with the main ion flow that it usually competes with, and the balance of
the two has a lot to say about whether impurities stay near a target or
leak upstream.

Flan does not add this force as a term in the equation of motion. It comes
out of the [collision model](collisions.md), once the background particles
an impurity collides with are sampled from the right distribution. This
page explains how, where the needed heat fluxes come from, and how well
it works.

!!! note "Where it is available"
    The thermal force needs the parallel heat fluxes of the background.
    Flan currently gets them from SOLEDGE3X backgrounds only, and with them
    the collisions run on the CPU only. For other backgrounds a warning at
    the start of the run says the thermal force isn't included.

## Where the force comes from

Coulomb collisions are much stronger with slow particles than with fast
ones: the collision rate goes as $1/g^3$, where $g$ is the relative speed.

In a plasma with a temperature gradient, the particles arriving at a point
from the hot side are faster than those arriving from the cold side. Since
there is no net flow (beyond the plasma's flow velocity), the slower
particles from the cold side must be slightly more numerous to make up for
it. An impurity sitting there is therefore hit from both sides, but the
particles from the cold side are slower and collide with it more
effectively. The net push is toward the hot side.

The same lopsidedness of the distribution is what carries heat down the
temperature gradient. So the thermal force and the conductive heat flux
are two faces of one thing: the background distribution isn't a
Maxwellian. A Maxwellian is symmetric about the flow velocity, carries no
heat flux, and gives no thermal force, no matter how the temperature
varies in space.

This is why the Nanbu collision model, which samples the colliding
background particle from a Maxwellian, has no thermal force by itself, and
also how to give it one: sample from the skewed distribution instead.

## The skewed distribution

With `imp_thermal_force = "on"` (the default), the background ion or
electron that an impurity collides with is sampled from

$$
f = f_M \left[1 + \epsilon\, w_\parallel \left(a + b\,w^2 + c\,w^4\right)\right],
$$

where $f_M$ is the Maxwellian with the local density, flow and
temperature, $\mathbf{w} = (\mathbf{v} - \mathbf{u})/v_\text{th}$ is the
velocity relative to the flow in units of the thermal speed
$v_\text{th} = \sqrt{T/m}$, and $w_\parallel$ is its component along the
magnetic field. Two things describe the skew.

### Its size, $\epsilon$

$$
\epsilon = \frac{q_\parallel}{n\,T\,v_\text{th}},
$$

where $q_\parallel$ is the species' parallel conductive heat flux,
positive along $\mathbf{B}$. $n T v_\text{th}$ is the plasma's thermal
energy density moving at the thermal speed, roughly the largest heat flux
it could carry, so $\epsilon$ is the heat flux as a fraction of that.
$\epsilon = 0$ is a Maxwellian, and $|\epsilon| \sim 0.1$ is a mild skew.
Each species has its own $\epsilon$.

### Its shape, $a$, $b$ and $c$

The shape says at which speeds the particles lean along the heat flux and
at which against it. It is set by three conditions on the skewed
distribution:

1. **No extra flow.** The flow is already in the Maxwellian's drift, so
   the skew must not add any.
2. **The heat flux is $q_\parallel$.**
3. **The push on a slow, heavy impurity is right.** The push is
   proportional to the average of $w_\parallel^2 (a + b w^2 + c w^4)/w^3$
   over the Maxwellian, which weights slow particles heavily.

The first two alone give the lowest-order shape, $w_\parallel (w^2/5 - 1)$
(Grad's 13-moment distribution, $c = 0$). It carries the right heat flux
but puts too much of the lean on slow particles, so it pushes too hard.
With the heat flux of a collisional plasma, $q_\parallel =
-\kappa \nabla_\parallel T$ with Braginskii's conductivities, it would
give $3.31\,Z^2 \nabla_\parallel T_i$ from the ions and
$1.90\,Z^2 \nabla_\parallel T_e$ from the electrons.

The third condition fixes this. Writing the wanted push as a fraction
"push" of the lowest-order shape's, the three conditions give

$$
c = \frac{1 - \text{push}}{25}, \qquad b = \frac{1}{5} - 14c, \qquad a = -1 + 35c .
$$

| Species | Push | Thermal force it gives |
|---|---|---|
| Ions | 0.80 | $2.65\,Z^2\,\nabla_\parallel T_i$ |
| Electrons | 0.37 | $0.71\,Z^2\,\nabla_\parallel T_e$ |

A push below 1 moves the lean to faster particles, which collide less. The
electron value is much lower because collisions with ions slow fast
electrons far less than slow ones, so the real distortion of the electron
distribution sits at high speeds.

The push is a property of the background's distribution, not of the
impurity. The values above are those of a collisional (fluid) plasma, and
the electron one is for a hydrogenic plasma ($Z_\text{eff} = 1$).

## Where the heat fluxes come from

For a SOLEDGE3X background, Flan computes $q_\parallel$ for the ions and
electrons from the fluxes SOLEDGE3X saves on the cell faces (`fluxn` and
`fluxE` in the plasma files), which are the fluxes its own equations use.
No temperature gradient is taken and no conductivity is assumed.

- **Only the conductive part is used.** `fluxE` is the total energy flux.
  The part carried by the flow, $(\tfrac{5}{2}T + \tfrac{1}{2}m u_\parallel^2)$
  times the particle flux, is subtracted, because the flow is already in
  the Maxwellian's drift. What is left is the heat flux that goes with the
  skew.
- **The fluxes on the toroidal ($\phi$) faces are used.** The magnetic
  field is mostly toroidal, so the parallel fluxes cross these faces,
  while the fluxes on the poloidal faces are mixed with the cross-field
  transport. The flux through a $\phi$ face is the parallel flux times
  $b_\phi$ and the face's area, which is divided out. A cell's value is
  the average of its two $\phi$ faces.
- **The normalization is checked.** One constant, $C$, converts the saved
  fluxes to physical units. Flan fits it by comparing the saved particle
  flux with the parallel flow it already reads, prints it at the start of
  the run for each species, and compares it with the expected
  $-2\pi c_0/a_0$. A difference of more than 5% stops the run, since the
  heat fluxes couldn't be trusted.

The plasma files must have the face fluxes for both species, in every
frame. If no frame has them, the run goes ahead without the thermal force
and says so. If only some do, the run stops.

$q_\parallel$ of each species is saved with the background in the output
file (`q_par_i` and `q_par_e`, in W/m²).

## Sampling

A velocity is drawn from the skewed distribution by rejection: draw from
the Maxwellian, and keep the draw with probability $F/M$, where
$F = f/f_M$ is the skew factor at that velocity and
$M = 1 + |\epsilon| \times \text{bound}$ is an upper limit on it for speeds
up to $5\,v_\text{th}$. On average it takes $M$ draws per sample. Where
$q_\parallel = 0$, a single Maxwellian draw is made, exactly as with the
option off.

The expansion behind the skewed distribution assumes $\epsilon$ is small,
and three things happen where it isn't. All are counted.

| What | When | What is done |
|---|---|---|
| Clipped draw | $F < 0$: the skew factor is negative in the tails | The draw is dropped |
| Above the bound | $F > M$: a draw faster than $5\,v_\text{th}$ | The draw is kept |
| Capped skew | $\lvert\epsilon\rvert > 1$ | $\epsilon$ is set to $\pm 1$ |

Clipping makes the sampled distribution carry a little less heat flux
than asked for. For the ion shape, the sampled heat flux is 99.6% of the
target at $\epsilon = 0.07$ and 86% at $\epsilon = 0.3$, where 0.3% of the
draws are clipped. For the electron shape it is 94% at $\epsilon = 0.1$.
The density, flow and temperature are hardly affected (the temperature
rises by 0.7% at $\epsilon = 0.3$).

The counts are printed at the end of the run and saved in the
`thermal_force` group of the output file, for ions and electrons:

| Attribute | Meaning |
|---|---|
| `samples` | Background particles sampled with a skew |
| `draws` | Maxwellian draws, kept or not |
| `negative` | Clipped draws |
| `above_bound` | Draws above the bound |
| `eps_capped` | Samples where $\lvert\epsilon\rvert$ was capped at 1 |

`negative / draws` is the share of the Maxwellian that was clipped. If it
is more than a percent or so, or if many samples were capped, the heat
flux is too large for the expansion in much of the region the impurities
are in, and the thermal force there is only approximate.

## Tests

The tests in `tests/test_thermal_force_cpu.cpp` use the slab background
(`bkg_source = "slab"`): a uniform plasma with the magnetic field along
$x$ and ion and electron heat fluxes along the field with
$\epsilon = 0.1$ (`slab_heat_flux_i` and `slab_heat_flux_e`). The
regression case `regression/reg_thermal_force.cpp` checks the drift at
force balance in a full run.

- **The kick is the right size.** Tungsten at rest is collided millions
  of times, and its mean change in velocity along the field is compared
  with the exact expectation for Nanbu's model at that time step. It
  agrees within the statistical error (2 to 3%) for ions and electrons, at
  a time step long enough ($s_0 = 3$) that the kick is well below its
  small-step value. The kick is against the heat flux, there is none
  across the field, and none with the option off.
- **The small-step limit is the textbook force.** With Braginskii's heat
  flux, the force in the limit of a small time step is $2.65\,Z^2
  \nabla_\parallel T_i$ and $0.71\,Z^2 \nabla_\parallel T_e$, to 0.2%.
  This is a calculation, not a simulation: it ties together the collision
  rate in Nanbu's $s$, the conductivities and the fitted shapes.
- **Force balance.** Tungsten left to collide with ions that carry a heat
  flux and don't flow doesn't come to rest. It drifts against the heat
  flux at the speed where friction cancels the thermal force, and that
  speed matches the expected one within the statistical error (2%). Its
  temperature settles at the ions'.
- **Sub-cycling recovers the force.** One time step with $s_0 \approx 3$
  is split into 30 sub-steps, and the kick matches 30 times the exact
  kick of one such sub-step, well above the kick without sub-cycling. A
  step that doesn't need splitting gives exactly the same result as with
  sub-cycling off.

## The time step

Nanbu's model turns the relative velocity of each sampled pair by an
angle that depends on $s \propto \Delta t / g^3$ of that pair (see
[Collisions](collisions.md#the-nanbu-model)). The thermal force comes
from the slowest background particles, which have the largest $s$. Once
their $s$ is large they are already deflected as much as they can be, and
a longer time step adds nothing more. The force per unit time therefore
falls as the time step grows.

Let $s_0$ be $s$ for a background particle at the thermal speed hitting an
impurity at rest. For ions:

| $s_0$ | Thermal force, relative to a small step | Drift at force balance, in $\epsilon\,v_\text{th}$ |
|---|---|---|
| $\to 0$ | 1.00 | $-0.48$ |
| 0.1 | 0.83 | $-0.45$ |
| 1 | 0.51 | $-0.39$ |
| 3 | 0.30 | $-0.32$ |
| 10 | 0.11 | $-0.22$ |

The friction on a slow impurity falls too, though less, so the drift at
which the two balance is less affected than the force.

For the same time step, $s_0$ of the electrons is about
$\sqrt{m_i/m_e} \approx 60$ times that of the ions. In hot, thin plasma
both are small. In cold, dense plasma the ion value is usually still
small, but the electron one may not be, and the electron thermal force is
then underestimated. $s_0$ grows with density and $Z^2$, and falls with
temperature as $T^{-3/2}$.

### Sub-cycling

A smaller `imp_time_step` would fix this, but it slows down everything
else. Instead, Flan splits only the collision. When $s_0$ of a time step
is above a target, the collision is done as $n$ collisions of
$\Delta t / n$, each with a newly sampled background particle:

$$
n = \min\left(n_\text{cap},\ \left\lceil \frac{s_0}{s_\text{target}} \right\rceil\right)
$$

The impurity then meets $n$ background particles in the step in place of
one, and each sub-step has $s_0 / n$. The plasma at the impurity's
position is interpolated once for all of them, and the impurity's
velocity is updated after each. Ions and electrons are handled
separately, each with its own $s_0$. Where $s_0$ is already at or below
the target, nothing changes.

The target sets how close the forces get to their small-step values, as
long as the cap isn't reached:

| Target | Ion thermal force | Electron thermal force | Friction |
|---|---|---|---|
| 0.1 | 0.83 | 0.91 | 0.88 |
| 0.03 (default) | 0.91 | 0.96 | 0.94 |
| 0.01 | 0.95 | 0.98 | 0.97 |

The shortfall shrinks slowly, as about $s_\text{target}^{2/3}$, because
there are always some background particles slow enough to have $s > 1$
(those below $s_0^{1/3}$ thermal speeds), and the forces lean heavily on
the slow ones. Getting within a percent would take a target near 0.001.

The cap limits the cost in very collisional plasma. Once $s_0$ exceeds
$s_\text{target}\,n_\text{cap}$ (3.84 with the defaults), the forces fall
off again as in the table above, with $s_0 / n_\text{cap}$ in place of
$s_0$. An impurity in plasma that collisional is already tied closely to
the background flow.

The two options are `imp_coll_subcycle_target` (default 0.03) and
`imp_coll_subcycle_cap` (default 128). A cap of 1 turns sub-cycling off.
It applies to all Nanbu collisions, with or without the thermal force,
and on the CPU and the GPU.

## Using it

The thermal force is on by default, so with a SOLEDGE3X background whose
plasma files have the face fluxes, and `use_gpu = "off"`, nothing needs
to be set. To compare with a run without it:

```cpp
inpts["imp_thermal_force"] = "off";
```

What Flan does in each situation:

| Situation | Result |
|---|---|
| Heat fluxes available, CPU | Thermal force included |
| No heat fluxes (other backgrounds, or files without face fluxes) | Warning, run goes ahead without it |
| Heat fluxes available, GPU | Run stops: use the CPU or turn the option off |
| Only one species has a heat flux | Run stops |
| Option off, collisions on | Warning that the thermal force is left out |

## Limitations

- **Parallel only.** The skew is along the magnetic field. The thermal
  force from temperature gradients across the field is not included.
- **Average force only.** The skew is fitted to give the right mean force
  on a slow, heavy impurity. For lighter impurities, whose thermal speed
  isn't small next to the background ions', the collisions still include
  the impurity's own motion, but the result hasn't been checked against
  the mass-dependent coefficients in the literature.
- **Fitted to a collisional plasma.** The push values are those of a
  fluid closure. Where the plasma is not collisional the heat flux read
  from the background is still used as is, but the shape of the real
  distribution may differ. A kinetic background could supply both the
  heat flux and the push from its own distribution.
- **Hydrogenic electron coefficient.** The electron push is for
  $Z_\text{eff} = 1$.
- **Large skews are approximate.** See [Sampling](#sampling).
- **A few percent low.** Even with sub-cycling the force is somewhat
  below the fluid value, and more so where the cap on sub-steps is
  reached. See [The time step](#the-time-step).
- **SOLEDGE3X and CPU only**, for now.

## References

- S. I. Braginskii, "Transport processes in a plasma", in *Reviews of
  Plasma Physics*, Vol. 1, p. 205 (1965).
- S. Chapman, "Thermal diffusion in ionized gases", *Proc. Phys. Soc.*
  **72**, 353 (1958).
- H. Grad, "On the kinetic theory of rarefied gases", *Commun. Pure Appl.
  Math.* **2**, 331 (1949).
- K. Nanbu, "Theory of cumulative small-angle collisions in plasmas",
  *Phys. Rev. E* **55**, 4642 (1997).
- P. C. Stangeby, *The Plasma Boundary of Magnetic Fusion Devices*, IOP
  Publishing (2000).
