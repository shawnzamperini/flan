# Impurity Transport

Flan follows impurity ions through a background plasma (usually from Gkeyll) by 
integrating each particle's full orbit, including its gyromotion. No 
guiding-center approximation is made. Each particle is pushed through the 
background electric and magnetic fields with the Boris algorithm, and its
charge state and velocity are updated between pushes by ionization/recombination
and Coulomb collisions.

## The main loop

Particles are held in a fixed number of *slots*. Each pass through the main loop
(`ImpurityTransport` namespace, `impurity_transport.cpp`) advances every alive
particle in the slots by one time step, `imp_time_step`:

1. **Ionization/recombination** (if `imp_iz_recomb = "on"`): the charge state
   may change based on OpenADAS rate coefficients. See
   [Ionization/Recombination](ionization-recombination.md).
2. **Collisions** (if `imp_collisions = "nanbu"`): the velocity is scattered
   by the background ions and electrons with the Nanbu model. See
   [Collisions](collisions.md).
3. **Boris push**: the velocity is updated from the Lorentz force
   (described below).
4. **Step**: the position is moved forward by $\mathbf{v}\,\Delta t$, and the
   particle's time is advanced.
5. **Bounds check**: particles that leave the domain are absorbed, wrapped
   (periodic boundaries) or otherwise handled by the boundary conditions.
6. **Fill slots**: the slots of dead particles are filled with newly launched
   particles until none are left to follow.
7. **Find cell**: the grid cell and background time frame that contain each
   particle are found.
8. **Record statistics**: each particle's contribution to the impurity density
   and other outputs is added to the cell it's in.

The loop runs until every slot is dead and there are no particles left to 
launch. A particle is also killed once it has been followed for longer than
`imp_max_time`.

Each step has a CPU (OpenMP) and GPU (CUDA) implementation. A wrapper function
picks between them based on `use_gpu`, so the physics is the same either way.

## Coordinates

The background is defined on a computational grid $(x, y, z)$ that is generally 
curvilinear, e.g., field-aligned coordinates. The Lorentz force is much simpler to 
evaluate in Cartesian space, though, so each particle carries two velocities:

- $(v_X, v_Y, v_Z)$: the Cartesian velocity. This is what the Boris algorithm
  updates.
- $(v_x, v_y, v_z)$: the velocity in computational coordinates. This is what moves
  the particle on the grid.

After the Cartesian velocity has been updated, the computational velocity is found with the
reciprocal (contravariant) basis vectors $\nabla x$, $\nabla y$, $\nabla z$:

$$
v_x = \frac{\partial x}{\partial X} v_X + \frac{\partial x}{\partial Y} v_Y 
    + \frac{\partial x}{\partial Z} v_Z
$$

and similarly for $v_y$ and $v_z$. The basis vectors are stored at the cell 
centers, and are trilinearly interpolated to the particle's position.

## The Boris algorithm

The Boris algorithm is the standard method for pushing charged particles
through electric and magnetic fields. It is second-order accurate and 
time-reversible. A pure magnetic field does no work, and the algorithm
conserves energy exactly in that case, so the
gyro-orbit doesn't drift over the many gyro-periods a particle is followed for.

### Leapfrog staggering

Positions and velocities are staggered by half a time step (leapfrog):

$$
\mathbf{x}_n = \mathbf{x}(t_n), \qquad \mathbf{v}_n = \mathbf{v}(t_n - \Delta t/2)
$$

Particles are launched with the position and velocity at the same time, so the
velocity is pushed *backwards* by half a time step ($-\Delta t/2$) before the
main loop starts. This sets up the staggering.

### Field interpolation

Before the push, $\mathbf{E}$ and $\mathbf{B}$ are found at the particle's position and time:

1. Each Cartesian component, and the magnitude, is trilinearly interpolated in
   space between the centers of the cell containing the particle and its
   neighbors on the particle's side. This is done for the particle's
   background time frame and the next one (the previous one if the particle is in the last
   frame).
2. The two values are linearly interpolated in time to the particle's time.
3. The interpolated vector is rescaled to the interpolated magnitude.

Interpolation keeps the fields continuous as a particle crosses cell boundaries
and background time frames. With piecewise-constant fields, particles get
unphysical "kicks" at every boundary and drift off the field line (numerical
diffusion). Interpolating in time also matters for time-dependent drifts like the 
polarization drift.

The rescaling in step 3 is a deliberate trade-off. An interpolated vector is
almost always shorter than the true field, which would make the gyroradius
wrong. Rescaling to the interpolated magnitude gets the gyro-orbit right, at the cost of a
small discontinuity at cell boundaries. A particle spends far more time inside
a cell than crossing a boundary, so getting the orbit right matters more.

### The push

With the charge-to-mass ratio $q/m$ (where $q = Z e$ for charge state $Z$), the 
velocity is advanced from $\mathbf{v}_n$ to $\mathbf{v}_{n+1}$:

1.  **Half electric kick:**

    $$
    \mathbf{v}^- = \mathbf{v}_n + \frac{q}{m}\mathbf{E}\,\frac{\Delta t}{2}
    $$

2.  **Magnetic rotation.** With

    $$
    \mathbf{t} = \frac{q}{m}\mathbf{B}\,\frac{\Delta t}{2}, \qquad
    \mathbf{s} = \frac{2\mathbf{t}}{1 + |\mathbf{t}|^2}
    $$

    the velocity is rotated about $\mathbf{B}$ in two steps:

    $$
    \mathbf{v}' = \mathbf{v}^- + \mathbf{v}^- \times \mathbf{t}, \qquad
    \mathbf{v}^+ = \mathbf{v}^- + \mathbf{v}' \times \mathbf{s}
    $$

3.  **Second half electric kick:**

    $$
    \mathbf{v}_{n+1} = \mathbf{v}^+ + \frac{q}{m}\mathbf{E}\,\frac{\Delta t}{2}
    $$

The rotation leaves $|\mathbf{v}^+| = |\mathbf{v}^-|$, so all energy change comes
from the electric field. Finally, $(v_x, v_y, v_z)$ are calculated from the new
Cartesian velocity (see [Coordinates](#coordinates)), and the step moves the
particle:

$$
\mathbf{x}_{n+1} = \mathbf{x}_n + \mathbf{v}_{n+1}\,\Delta t
$$

### Implementation

The CPU version is `Boris::update_velocity_cpu` in `src/boris.cpp`, parallelized
over particles with OpenMP. The GPU version is `update_velocity_kernel` in
`src/cuda/boris.cu`, with one thread per particle. The two use the same push.
The GPU version builds the interpolation stencils with the helpers in
`src/cuda/interpolate.cu`. The core of the push is:

```cpp
// Charge state times the (positive) elementary charge, divided by mass
double q_m {slots.q()[i] * -Constants::charge_e / slots.mass()};

// t vector
for (int j {}; j < 3; ++j) t[j] = q_m * B_local[j] * 0.5 * dt;
double tmag2 {t[0]*t[0] + t[1]*t[1] + t[2]*t[2]};

// s vector
for (int j {}; j < 3; ++j) s[j] = 2 * t[j] / (1.0 + tmag2);

// Half electric kick
for (int j {}; j < 3; ++j) vminus[j] = v[j] + q_m * E_local[j] * 0.5 * dt;

// Magnetic rotation
vminus_cross_t = Utilities::cross_product(vminus, t);
for (int j {}; j < 3; ++j) vprime[j] = vminus[j] + vminus_cross_t[j];
vprime_cross_s = Utilities::cross_product(vprime, s);
for (int j {}; j < 3; ++j) vplus[j] = vminus[j] + vprime_cross_s[j];

// Second half electric kick
slots.set_vX(i, vplus[0] + q_m * E_local[0] * 0.5 * dt);
slots.set_vY(i, vplus[1] + q_m * E_local[1] * 0.5 * dt);
slots.set_vZ(i, vplus[2] + q_m * E_local[2] * 0.5 * dt);
```

`Constants::charge_e` is the electron charge ($-1.602\times10^{-19}$ C), so
`q * -charge_e` is the positive ion charge $Ze$.

## Neutrals

!!! note "Neutrals are handled correctly"
    Particles in charge state $Z = 0$ are followed in the same loop as the ions and
    need no special treatment in the Boris push. With $q = 0$:

    - $q/m = 0$, so both electric kicks vanish,
    - $\mathbf{t} = \mathbf{s} = 0$, so there is no magnetic rotation,

    and $\mathbf{v}_{n+1} = \mathbf{v}_n$. A neutral moves in a straight line
    at constant velocity, as it should, and it isn't affected by the fields.

The other steps also treat neutrals correctly:

- **Collisions:** neutrals are skipped, since they don't undergo Coulomb collisions.
- **Ionization/recombination:** a neutral can ionize to $Z = 1$ but can't
  recombine, so its recombination probability is zero. Once it's ionized, the
  next Boris push uses the new charge, and the particle starts to gyrate.

A neutral's velocity can only change when it ionizes. Until then, it 
just streams through the plasma.
