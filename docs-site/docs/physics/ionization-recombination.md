# Ionization and Recombination

As an impurity moves through the plasma, its charge state changes through
electron-impact ionization and recombination. Flan models both as random events
whose probabilities come from OpenADAS rate coefficients. They are turned on with
`imp_iz_recomb = "on"`, and are applied once per time step, at the start of
each pass through the main loop (see
[Impurity Transport](impurity-transport.md#the-main-loop)).

## Rate coefficients

The rate coefficients come from the OpenADAS database
([open.adas.ac.uk](https://open.adas.ac.uk)), in ADF11 format:

- **SCD** files: effective ionization rate coefficients, $S_Z(n_e, T_e)$, for
  $Z \to Z+1$.
- **ACD** files: effective recombination rate coefficients,
  $\alpha_{Z}(n_e, T_e)$, for $Z \to Z-1$.

These are *effective* (collisional-radiative) coefficients. They depend on
both $n_e$ and $T_e$, and they include the effect of the excited-state populations.
Recombination includes radiative, dielectronic and three-body
recombination.

Each file is read from

```
[openadas_root]/scd[openadas_year]/scd[openadas_year]_[element].dat
[openadas_root]/acd[openadas_year]/acd[openadas_year]_[element].dat
```

where `openadas_root` and `openadas_year` are input options (`openadas_year`
defaults to `50`) and the element is chosen with `imp_atom_num`. For example,
tungsten with the default year reads `scd50/scd50_w.dat` and `acd50/acd50_w.dat`.
The supported elements are H through Ne, Al, Si, Cl, Ar, Fe, Kr, Mo, Xe and W.

The files store $\log_{10}$ of the densities (cm$^{-3}$), temperatures (eV) and
rate coefficients (cm$^3$/s). These are converted to m$^{-3}$, eV and m$^3$/s
when they are read in (`OpenADAS::read_rate_coefficients` in `src/openadas.cpp`).

### Interpolation

For each particle, $n_e$ and $T_e$ are interpolated to the particle's position and
time, the same way as in the [Boris push](impurity-transport.md#field-interpolation).
The rate coefficient for the particle's charge state is then bilinearly
interpolated in $(T_e, n_e)$ from the table.

If $n_e$ or $T_e$ fall outside the range of the OpenADAS data, they are clamped
to the nearest edge of the table. The rates are never extrapolated.

## Charge-state update

During one time step $\Delta t$, a particle in charge state $Z$ ionizes with
probability

$$
P_\text{iz} = S_Z(n_e, T_e)\, n_e\, \Delta t
$$

and recombines with probability

$$
P_\text{rec} = \alpha_Z(n_e, T_e)\, n_e\, \Delta t.
$$

Two uniform random numbers $U_1, U_2 \in [0, 1)$ are drawn, and the charge
becomes

$$
Z \to Z + [U_1 < P_\text{iz}] - [U_2 < P_\text{rec}],
$$

where $[\cdot]$ is 1 if the condition is true and 0 otherwise. If both events
happen in the same step they cancel and the charge state is unchanged, so the
charge changes by at most one per step.

A fully stripped ion ($Z$ equal to the atomic number) can't ionize, so
$P_\text{iz} = 0$. The velocity is not changed by an ionization or recombination
event. The new charge simply takes effect in the next Boris push and collision
step.

The CPU version is `OpenADAS::ioniz_recomb` (`src/openadas.cpp`), and the GPU
version is `ioniz_recomb_kernel` (`src/cuda/openadas.cu`).

### Neutrals

!!! note "Neutrals are handled correctly"
    A neutral ($Z = 0$) can ionize to $Z = 1$ using the $Z = 0 \to 1$ SCD
    coefficient, but it can't recombine, so $P_\text{rec} = 0$. Once it's
    ionized, it is pushed by the fields and starts to undergo collisions. See
    [Impurity Transport](impurity-transport.md#neutrals) and
    [Collisions](collisions.md#neutrals).

## Choosing the time step

$P = R\, n_e\, \Delta t$ is the first-order expansion of the exact
probability $1 - e^{-R\, n_e\, \Delta t}$, and it's only accurate when $P \ll 1$.
The update also allows at most one ionization or recombination per step. Both
mean the time step must be short compared to the ionization and recombination
times, $1/(S_Z n_e)$ and $1/(\alpha_Z n_e)$.

If $P > 1$, the event always happens and the true rate is underestimated. Flan
counts how often this happens, and prints a warning at the end of the run:

```
Warning: Ionization probability was > 1.0 for N particles. Consider using a smaller time step.
```

(and similarly for recombination). The counter is incremented each time step that the
probability is above 1, so $N$ counts particle-steps rather than unique particles. If these warnings
show up for more than a handful of particle-steps, reduce `imp_time_step`. This is most
likely for low charge states in hot, dense plasma, where ionization is fastest.

## When the model is valid

**Valid when:**

- **The time step resolves the atomic time scales**, i.e., $P \ll 1$ (see above).
- **The plasma is within the range of the OpenADAS tables.** Outside the table,
  the edge values are used, which can be a poor approximation, e.g., at
  temperatures far below the lowest tabulated $T_e$.
- **The electrons are close to Maxwellian.** The ADAS coefficients are
  averaged over a Maxwellian electron distribution at $T_e$.
- **The excited states are in quasi-steady state.** The effective coefficients
  assume that excited-state populations relax much faster than the plasma
  conditions change along the impurity's path. This usually holds, though
  long-lived metastable states can violate it.

**Not included:**

- **Charge-exchange recombination** with background neutrals, e.g., $W^{Z+} + D^0
  \to W^{(Z-1)+} + D^+$, which can be important in regions with a high neutral density.
- **Ionization by ions or photons.** Only electron-impact ionization is
  included.
- **More than one charge change per step.** See
  [Choosing the time step](#choosing-the-time-step).
