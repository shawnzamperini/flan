#!/usr/bin/env python3
"""
Checks the collision regression cases against fluid (analytic) results and
plots them.

Cases (each made by the program of the same name in regression/):

  reg_friction_<variant>.nc
      Tungsten starting at rest in a plasma flowing along B. The fluid
      result is <v_par>(t) = u (1 - exp(-t / tau_s)), with u the flow
      speed and tau_s the Spitzer slowing-down time.

  reg_thermal_force_<variant>.nc
      Tungsten starting at rest in a plasma at rest that carries parallel
      heat fluxes. The fluid result is a drift against the heat flux at
      the speed where friction cancels the thermal force,
      v_drift = -Z^2 (beta_i grad T_i + 0.71 grad T_e) tau_s / m_z,
      with grad T = q_par / kappa (Braginskii) and Stangeby's
      mass-dependent beta_i (2.65 for a heavy impurity).

  Both also have a light-impurity variant (C4+ starting at the ion
  temperature). Its friction rate is that between two Maxwellians, and
  the collision model's expectation (check 2) is only worked out for a
  slow impurity, so a light one is checked against fluid only.

Every number the checks need (density, temperatures, flow, masses, charge,
time step, sub-cycling settings) is read from the inputs saved in the
output file, so what is checked is always what was run.

Each variant is checked twice:

  1. Against the fluid result, to within TOL_FLUID (10%). This is the
     physics check. A miss is reported but not counted as a failure when
     most of the force being checked comes from a species whose
     collisions the run doesn't resolve: its sub-steps are still more
     collisional than imp_coll_subcycle_target, because sub-cycling is
     off or the cap on sub-steps (imp_coll_subcycle_cap) is reached.
     Those variants exist to show what sub-cycling corrects, or what the
     cap costs, and are expected to miss.
  2. Against what Flan's collision model (Nanbu's) should give at the time
     step and sub-cycling of the run, to within TOL_MODEL (3%) plus the
     statistical error. The model is a little below the fluid result at
     any finite step (see the thermal force page of the docs), by an
     amount that can be calculated. This tighter check catches a change
     in the code that the 10% check would let through.

After the checks, the effective thermal-force coefficients (beta_i,
alpha_e) that the thermal-force runs correspond to are printed for
comparison with the fluid and Stangeby values (see thermal_coefficients).
They are not checked.

Usage:
    python check_collisions.py [directory with the .nc files]

Writes check_friction.png and check_thermal_force.png to that directory,
prints a table, and exits
with status 1 if a check failed or no output files were found.
"""

import glob
import math
import os
import sys

import numpy as np
import netCDF4

# ----------------------------------------------------------------------
# Constants and tolerances
# ----------------------------------------------------------------------
E_CHARGE = 1.602176634e-19   # C
EPS0 = 8.8541878128e-12      # F/m
AMU = 1.66053906660e-27      # kg

TOL_FLUID = 0.10   # agreement with the fluid result
TOL_MODEL = 0.03   # agreement with the collision model's own expectation
N_SIGMA = 3.0      # statistical errors allowed on top of a tolerance


# ----------------------------------------------------------------------
# Reading a Flan output file
# ----------------------------------------------------------------------
class Run:
    """Inputs and the mean impurity velocity against time of one run."""

    def __init__(self, path):
        self.path = path
        self.name = os.path.splitext(os.path.basename(path))[0]
        with netCDF4.Dataset(path) as nc:
            inp = nc["input"]

            def scalar(key):
                return float(np.ravel(inp[key][:])[0])

            self.ne = scalar("slab_ne")                    # m^-3
            self.te = scalar("slab_te")                    # eV
            self.ti = scalar("slab_ti")                    # eV
            self.flow = scalar("slab_flow")                # m/s, along +X
            self.eps_i = scalar("slab_heat_flux_i")
            self.eps_e = scalar("slab_heat_flux_e")
            self.m_i = scalar("gkyl_ion_mass_amu") * AMU   # kg
            self.m_e = scalar("gkyl_elec_mass_amu") * AMU  # kg
            self.m_imp = scalar("imp_mass_amu") * AMU      # kg
            self.charge = int(scalar("imp_init_charge"))
            # Starting temperature (eV): 0 for impurities starting at
            # rest (the heavy cases), the ion temperature for a light
            # impurity, which would otherwise heat up during the run
            self.t_imp = scalar("imp_temp_start_val")
            self.dt = scalar("imp_time_step")              # s
            self.imp_num = int(scalar("imp_num"))
            self.sub_target = scalar("imp_coll_subcycle_target")
            self.sub_cap = int(scalar("imp_coll_subcycle_cap"))
            self.thermal_force = "".join(
                np.ravel(inp["imp_thermal_force"][:]).astype(str)
            ).strip().strip("\x00") == "on"

            # Mean velocity of all the impurities in each time frame: the
            # cell averages weighted by the counts in each cell. B is
            # along X, so v_X is the parallel velocity.
            self.t = np.array(nc["geometry"]["time"][:], dtype=float)
            counts = np.array(nc["output"]["Nz"][:], dtype=float)
            total = counts.sum(axis=(1, 2, 3))
            total = np.where(total > 0, total, np.nan)

            def mean(key):
                v = np.array(nc["output"][key][:], dtype=float)
                return (v * counts).sum(axis=(1, 2, 3)) / total

            self.v_par = mean("v_X")
            self.v_perp1 = mean("v_Y")
            self.v_perp2 = mean("v_Z")

    @property
    def subcycled(self):
        return self.sub_cap > 1

    @property
    def slow(self):
        """
        Whether the impurities are slow compared with the background ions:
        their thermal speed squared is under 5% of the ions'. The
        collision model's expectations (model_friction_rate, model_drift)
        are worked out for a slow impurity and are only checked then.
        """
        return self.m_i * self.t_imp / (self.m_imp * self.ti) < 0.05

    def v_err(self):
        """
        Statistical error of the mean velocity in a frame (m/s): the
        thermal spread of the impurities once they are at the ion
        temperature, over sqrt(number of particles). Before that the
        spread is smaller, so this is on the safe side.
        """
        return math.sqrt(self.ti * E_CHARGE / self.m_imp / self.imp_num)


# ----------------------------------------------------------------------
# Physics
# ----------------------------------------------------------------------
def coulomb_log(run):
    """
    Coulomb logarithm as Flan defines it: ln(Debye length / b0), with b0 =
    Z e^2 / (2 pi eps0 3 T) the mean impact parameter for a 90 degree
    deflection, and at least 1. The fluid formulas take the Coulomb
    logarithm as given, so the same one is used for them.
    """
    te_j = run.te * E_CHARGE
    debye = math.sqrt(EPS0 * te_j / (run.ne * E_CHARGE**2))

    def lnl(t_ev):
        b0 = run.charge * E_CHARGE**2 / (2.0 * math.pi * EPS0 * 3.0
                                         * t_ev * E_CHARGE)
        return max(math.log(debye / b0), 1.0)

    return lnl(run.ti), lnl(run.te)   # for ion and electron collisions


def spitzer_rate(run, m_b, t_b, lnl):
    """
    Slowing-down rate (1/s) of a test particle moving slowly through
    Maxwellian field particles of mass m_b and temperature t_b (eV), from
    the standard test-particle formulas (Spitzer; the slow limit of the
    NRL Plasma Formulary's nu_s):

        nu_s = (1 + m_z / m_b) (4 / (3 sqrt(pi))) (m_b / (2 T))^(3/2)
               n Z^2 e^4 lnL / (4 pi eps0^2 m_z^2)

    For impurities with a Maxwellian spread at temperature T_z, the rate
    at which their mean velocity relaxes is that between two Maxwellians:
    (m_b / 2 T)^(3/2) becomes (2 T / m_b + 2 T_z / m_z)^(-3/2), i.e. the
    rate above times (1 + m_b T_z / (m_z T))^(-3/2). This matters for a
    light impurity (about 0.79 for C in D at T_z = T_i); for tungsten
    starting at rest T_z = 0.
    """
    t_j = t_b * E_CHARGE
    spread = (1.0 + m_b * run.t_imp / (run.m_imp * t_b))**-1.5
    return ((1.0 + run.m_imp / m_b) * 4.0 / (3.0 * math.sqrt(math.pi))
            * (m_b / (2.0 * t_j))**1.5 * run.ne * run.charge**2
            * E_CHARGE**4 * lnl / (4.0 * math.pi * EPS0**2 * run.m_imp**2)
            * spread)


def fluid_friction_rate(run):
    """1 / tau_s: slowing down on the ions plus on the electrons."""
    lnl_i, lnl_e = coulomb_log(run)
    return (spitzer_rate(run, run.m_i, run.ti, lnl_i)
            + spitzer_rate(run, run.m_e, run.te, lnl_e))


def s0_of_step(run, m_b, t_b, lnl):
    """
    Nanbu's s for a background particle at its thermal speed, for one
    full time step: how collisional the step is (see collisions.h).
    """
    mu = run.m_imp * m_b / (run.m_imp + m_b)
    v_th = math.sqrt(t_b * E_CHARGE / m_b)
    sq = run.charge * E_CHARGE**2 / (EPS0 * mu)
    return lnl / (4.0 * math.pi) * sq * sq * run.ne * run.dt / v_th**3


def num_substeps(run, s0):
    """Sub-steps Flan splits a collision of this s0 into."""
    if not s0 > run.sub_target:
        return 1
    return int(min(run.sub_cap, math.ceil(s0 / run.sub_target)))


# Speeds (in thermal speeds) for averages over a Maxwellian
_W = np.logspace(-5.0, math.log10(9.0), 6000)
_FM = math.sqrt(2.0 / math.pi) * _W**2 * np.exp(-0.5 * _W**2)


def unresolved_share(run, weight_i, weight_e):
    """
    Share of a force that comes from a species whose collisions the run
    doesn't resolve, i.e., whose sub-steps have s0 above the target
    (sub-cycling off, or the cap reached). weight_i and weight_e are the
    ion and electron parts of the fluid force. The fluid check is only
    enforced when this is under a half.
    """
    lnl_i, lnl_e = coulomb_log(run)
    total = abs(weight_i) + abs(weight_e)
    if total == 0.0:
        return 0.0
    share = 0.0
    for m_b, t_b, lnl, weight in ((run.m_i, run.ti, lnl_i, weight_i),
                                  (run.m_e, run.te, lnl_e, weight_e)):
        s0 = s0_of_step(run, m_b, t_b, lnl)
        if s0 / num_substeps(run, s0) > 1.001 * run.sub_target:
            share += abs(weight) / total
    return share


def friction_factor(s0):
    """
    Mean fraction of a slow, heavy impurity's velocity (relative to the
    flow) that one Nanbu collision of this s0 removes, per unit of
    frac = m_b / (m_imp + m_b). A background particle at w thermal speeds
    has s = s0 / w^3 and turns the relative velocity so that <cos chi> =
    exp(-s); averaging over the Maxwellian gives <k - x exp(-x)>, with
    x = s0 / w^3 and k = 1 - exp(-x). For small s0 it is
    (1/3) sqrt(2 / pi) s0, which is the fluid friction.
    """
    x = s0 / _W**3
    return np.trapz((-np.expm1(-x) - x * np.exp(-x)) * _FM, _W)


def model_friction_rate(run):
    """
    Slowing-down rate (1/s) Flan's collision model should give at the
    time step and sub-cycling of the run: each sub-step removes
    frac * friction_factor(s0 / n) of the relative velocity.
    """
    lnl_i, lnl_e = coulomb_log(run)
    rate = 0.0
    for m_b, t_b, lnl in ((run.m_i, run.ti, lnl_i),
                          (run.m_e, run.te, lnl_e)):
        s0 = s0_of_step(run, m_b, t_b, lnl)
        n_sub = num_substeps(run, s0)
        frac = m_b / (run.m_imp + m_b)
        rate += frac * n_sub * friction_factor(s0 / n_sub) / run.dt
    return rate


# Thermal force. The background an impurity collides with is skewed by
# the heat flux: f = f_M [1 + eps w_par (a + b w^2 + c w^4)], eps = q_par /
# (n T v_th). The shape is set by its "push" (see collisions.h), fitted so
# that the force for small steps is the fluid closure's.
ALPHA_I, ALPHA_E = 2.6517, 0.71     # thermal force / (Z^2 grad T)
KAPPA_I, KAPPA_E = 3.906, 3.16      # conductivity / (n T tau / m)
PUSH_I = ALPHA_I / (0.6 * math.sqrt(2.0) * KAPPA_I)
PUSH_E = ALPHA_E / (0.6 * KAPPA_E)


def braginskii_time(run, elec, lnl):
    """Braginskii's collision time (s) of the ions or the electrons."""
    m_b, t_b = (run.m_e, run.te) if elec else (run.m_i, run.ti)
    return ((6.0 * math.sqrt(2.0) if elec else 12.0) * math.pi**1.5
            * EPS0**2 * math.sqrt(m_b) * (t_b * E_CHARGE)**1.5
            / (lnl * E_CHARGE**4 * run.ne))


def grad_t(run, elec):
    """
    Temperature gradient (J/m) that would carry the run's heat flux in a
    collisional plasma, |grad T| = q_par / kappa, with q_par = eps n T
    v_th and Braginskii's kappa. Its sign is eps's: the heat flows along
    +X, so the plasma is hotter toward -X.
    """
    lnl_i, lnl_e = coulomb_log(run)
    m_b, t_b, eps, lnl, coeff = (
        (run.m_e, run.te, run.eps_e, lnl_e, KAPPA_E) if elec
        else (run.m_i, run.ti, run.eps_i, lnl_i, KAPPA_I))
    t_j = t_b * E_CHARGE
    q_par = eps * run.ne * t_j * math.sqrt(t_j / m_b)
    kappa = coeff * run.ne * t_j * braginskii_time(run, elec, lnl) / m_b
    return q_par / kappa


def fluid_drift(run):
    """
    Drift (m/s) at which fluid friction cancels the fluid thermal force,
    -Z^2 (beta_i grad T_i + 0.71 grad T_e) / (m_z nu_s), with Stangeby's
    mass-dependent beta_i (stangeby_beta): 2.61 for W in D, 2.65 in the
    heavy limit, about 2.05 for C in D. Negative: against the heat flux.
    """
    if not run.thermal_force:
        return 0.0
    return -fluid_thermal_force(run) / (run.m_imp * fluid_friction_rate(run))


def fluid_thermal_force(run):
    """
    Fluid thermal force (N) along B, Z^2 (beta_i grad T_i + 0.71 grad
    T_e), whether or not the run has the thermal force on.
    """
    return run.charge**2 * (stangeby_beta(run) * grad_t(run, False)
                            + ALPHA_E * grad_t(run, True))


def thermal_factor(s0, push):
    """
    Mean kick along B that one Nanbu collision of this s0 gives a slow,
    heavy impurity, per unit of frac * v_th * eps: (1/3) <w^2 (a + b w^2 +
    c w^4) (1 - exp(-s0 / w^3))> over the Maxwellian. Negative: slow
    background particles, which lean against the heat flux, turn the
    furthest. For small s0 it is -0.2 sqrt(2 / pi) push s0.
    """
    c = (1.0 - push) / 25.0
    b = 0.2 - 14.0 * c
    a = -1.0 + 35.0 * c
    w2 = _W**2
    return np.trapz(w2 * (a + b * w2 + c * w2 * w2)
                    * -np.expm1(-s0 / _W**3) * _FM, _W) / 3.0


def model_drift(run):
    """
    Drift (m/s) Flan's collision model should give at the time step and
    sub-cycling of the run: where the thermal kicks and the friction of
    all the sub-steps of both species cancel.
    """
    if not run.thermal_force:
        return 0.0
    lnl_i, lnl_e = coulomb_log(run)
    kick, friction = 0.0, 0.0
    for m_b, t_b, lnl, eps, push in (
            (run.m_i, run.ti, lnl_i, run.eps_i, PUSH_I),
            (run.m_e, run.te, lnl_e, run.eps_e, PUSH_E)):
        s0 = s0_of_step(run, m_b, t_b, lnl)
        n_sub = num_substeps(run, s0)
        frac = m_b / (run.m_imp + m_b)
        v_th = math.sqrt(t_b * E_CHARGE / m_b)
        kick += frac * v_th * eps * n_sub * thermal_factor(s0 / n_sub, push)
        friction += frac * n_sub * friction_factor(s0 / n_sub)
    return kick / friction


# ----------------------------------------------------------------------
# Checks
# ----------------------------------------------------------------------
class Check:
    """One line of the results table."""

    def __init__(self, run, what, value, expected, tol, unit,
                 enforced=True):
        self.run, self.what, self.unit = run, what, unit
        self.value, self.expected, self.tol = value, expected, tol
        self.enforced = enforced
        self.ok = abs(value - expected) <= tol

    def status(self):
        if self.ok:
            return "pass"
        return "FAIL" if self.enforced else "miss (expected)"


def fit_rate(run):
    """
    Slowing-down rate (1/s) of the run: least-squares fit of
    u (1 - exp(-rate t)) to the mean parallel velocity, with u fixed at
    the flow speed. The first and last frames cover half a frame only
    and are left out. Returns the rate and its statistical error.
    """
    t, v = run.t[1:-1], run.v_par[1:-1]
    u = run.flow

    def residual(rate):
        return np.sum((v - u * -np.expm1(-rate * t))**2)

    # Golden-section search between generous bounds
    lo, hi = 0.01 / t[-1], 1000.0 / t[-1]
    lo, hi = math.log(lo), math.log(hi)
    gr = (math.sqrt(5.0) - 1.0) / 2.0
    for _ in range(200):
        a, b = hi - gr * (hi - lo), lo + gr * (hi - lo)
        if residual(math.exp(a)) < residual(math.exp(b)):
            hi = b
        else:
            lo = a
    rate = math.exp(0.5 * (lo + hi))

    # Error of the fitted rate from the error of each frame, v_err. The
    # frames are not independent: the mean velocity wanders over about a
    # slowing-down time, so frames closer than that move together. That
    # is allowed for by counting n_corr frames as one.
    sens = u * t * np.exp(-rate * t)                 # d(model) / d(rate)
    n_corr = max(1.0, 1.0 / (rate * (t[1] - t[0])))
    err = run.v_err() * math.sqrt(n_corr / np.sum(sens**2))
    return rate, err


def check_friction(run):
    checks = []
    rate, rate_err = fit_rate(run)
    fluid = fluid_friction_rate(run)
    model = model_friction_rate(run)
    run.rate, run.rate_fluid, run.rate_model = rate, fluid, model
    run.rate_err = rate_err

    # Slowing-down time, in microseconds, against the fluid result (not
    # enforced without sub-cycling) and against the model's expectation
    tau, tau_err = 1e6 / rate, 1e6 * rate_err / rate**2
    lnl_i, lnl_e = coulomb_log(run)
    resolved = unresolved_share(
        run, spitzer_rate(run, run.m_i, run.ti, lnl_i),
        spitzer_rate(run, run.m_e, run.te, lnl_e)) < 0.5
    checks.append(Check(run, "tau_s vs fluid", tau, 1e6 / fluid,
                        TOL_FLUID * 1e6 / fluid + N_SIGMA * tau_err, "us",
                        enforced=resolved))
    # The model's expectation is for a slow impurity: a light one is
    # checked against fluid only
    if run.slow:
        checks.append(Check(run, "tau_s vs model", tau, 1e6 / model,
                            TOL_MODEL * 1e6 / model + N_SIGMA * tau_err,
                            "us"))
    else:
        run.rate_model = float("nan")

    # Late in the run the impurities move with the flow. Over the last
    # quarter of the run, compare with u (1 - exp(-t / tau_s)) at the
    # model's rate (the fluid one for a light impurity), since a run that
    # slows down at the right rate should also end up in the right place
    # (the variant without sub-cycling hasn't fully caught up with the
    # flow by then).
    ref, tol = (model, TOL_MODEL) if run.slow else (fluid, TOL_FLUID)
    late = run.t[:-1] > 0.75 * run.t[-1]
    v_late = run.v_par[:-1][late].mean()
    expected = run.flow * (1.0 - np.exp(-ref * run.t[:-1][late]).mean())
    checks.append(Check(run, "final v_par", v_late, expected,
                        tol * abs(run.flow) + N_SIGMA * run.v_err(),
                        "m/s"))

    # Nothing across the field
    for label, v in (("mean v_Y", run.v_perp1), ("mean v_Z", run.v_perp2)):
        checks.append(Check(run, label, np.nanmean(v[1:-1]), 0.0,
                            N_SIGMA * run.v_err(), "m/s"))
    return checks


def check_thermal_force(run):
    checks = []
    fluid = fluid_drift(run)
    model = model_drift(run) if run.slow else float("nan")
    nu_s = fluid_friction_rate(run)
    run.drift_fluid, run.drift_model, run.rate_fluid = fluid, model, nu_s

    # Drift at force balance: the mean over frames later than 5 fluid
    # slowing-down times. The mean velocity wanders over about a
    # slowing-down time, so the frames average down as the number of
    # slowing-down times they span.
    late = run.t[:-1] > 5.0 / nu_s
    late[0] = False
    if late.sum() < 3:
        raise RuntimeError(f"{run.name}: run too short to reach force "
                           "balance")
    t_late = run.t[:-1][late]
    drift = run.v_par[:-1][late].mean()
    n_indep = max(1.0, (t_late[-1] - t_late[0]) * nu_s)
    drift_err = run.v_err() / math.sqrt(n_indep)
    run.drift, run.drift_err = drift, drift_err

    # The scale for the tolerances is the fluid drift with both heat
    # fluxes on, also for the variant with the thermal force off
    scale = abs(fluid_thermal_force(run) / (run.m_imp * nu_s))
    resolved = (not run.thermal_force) or unresolved_share(
        run, stangeby_beta(run) * grad_t(run, False),
        ALPHA_E * grad_t(run, True)) < 0.5
    checks.append(Check(run, "drift vs fluid", drift, fluid,
                        TOL_FLUID * scale + N_SIGMA * drift_err, "m/s",
                        enforced=resolved))
    if run.slow:
        checks.append(Check(run, "drift vs model", drift, model,
                            TOL_MODEL * scale + N_SIGMA * drift_err, "m/s"))

    # Nothing across the field
    for label, v in (("mean v_Y", run.v_perp1), ("mean v_Z", run.v_perp2)):
        checks.append(Check(run, label, np.nanmean(v[1:-1]), 0.0,
                            N_SIGMA * run.v_err(), "m/s"))
    return checks


# ----------------------------------------------------------------------
# Effective thermal-force coefficients
# ----------------------------------------------------------------------
def stangeby_beta(run):
    """
    Stangeby's mass-dependent ion thermal-force coefficient, divided by
    Z^2 so it compares with Flan's beta_i / Z^2:

        beta_i = 3 (mu + 5 sqrt(2) Z^2 (1.1 mu^(5/2) - 0.35 mu^(3/2)) - 1)
                 / (2.6 - 2 mu + 5.4 mu^2),   mu = m_z / (m_z + m_i).

    For a heavy impurity (mu -> 1, large Z) it goes to 2.65. It is given
    for reference: Flan's skew is fitted to the heavy limit.
    """
    mu = run.m_imp / (run.m_imp + run.m_i)
    z2 = run.charge**2
    beta = 3.0 * (mu + 5.0 * math.sqrt(2.0) * z2
                  * (1.1 * mu**2.5 - 0.35 * mu**1.5) - 1.0) \
        / (2.6 - 2.0 * mu + 5.4 * mu**2)
    return beta / z2


def same_plasma(a, b):
    """
    Whether two runs have the same plasma, impurity, time step and
    sub-cycling, so that the slowing-down rate measured in one (a
    reg_friction run) applies to the other. The flow doesn't matter: the
    friction depends on the velocity relative to it.
    """
    keys = ("ne", "te", "ti", "m_i", "m_e", "m_imp", "charge", "t_imp",
            "dt", "sub_target", "sub_cap")
    return all(math.isclose(getattr(a, k), getattr(b, k), rel_tol=1e-6,
                            abs_tol=1e-12) for k in keys)


class Coefficient:
    """One line of the effective-coefficient table."""

    def __init__(self, run, label, flan_fluid, flan, err, nu_from, model,
                 fitted, stangeby):
        self.run, self.label = run, label
        self.flan_fluid, self.flan, self.err = flan_fluid, flan, err
        self.nu_from, self.model = nu_from, model
        self.fitted, self.stangeby = fitted, stangeby


def thermal_coefficients(thermal, friction):
    """
    Effective thermal-force coefficients of the reg_thermal_force runs:
    the alpha_e and beta_i of F = Z^2 (beta_i grad T_i + alpha_e grad T_e)
    that Flan's drift at force balance corresponds to. At force balance
    the friction m_z nu_s v_drift cancels the thermal force, so

        coefficient = -m_z nu_s v_drift / (Z^2 grad T).

    A run with only the ion heat flux gives beta_i, one with only the
    electron heat flux alpha_e, and one with both (and equal temperature
    gradients, as in reg_thermal_force) their sum.

    The friction rate nu_s is taken two ways:

      - the fluid (Spitzer) rate, which gives the coefficient a fluid code
        would infer from Flan's drift. Nanbu's friction shortfall at a
        finite step shows up in it.
      - Flan's own rate: the one measured by a reg_friction run of the
        same plasma, time step and sub-cycling if there is one, else what
        the collision model should give (model_friction_rate; the fluid
        rate for an impurity that isn't slow). This isolates the thermal
        force itself.

    For comparison, the table also gives what the collision model should
    give at the run's step (its thermal kick per unit Z^2 grad T; slow
    impurities only), the heavy-limit value Flan's skew is fitted to, and
    Stangeby's mass-dependent beta_i, which is the fluid result the
    drift is checked against.

    grad T is not an input of the slab: it is the gradient that would
    carry the run's heat flux in a collisional plasma (grad_t), so these
    are coefficients per unit heat flux, converted with Braginskii's
    conductivities.
    """
    rows = []
    for run in thermal:
        if not run.thermal_force:
            continue
        g_i, g_e = grad_t(run, False), grad_t(run, True)
        if g_e == 0.0 and g_i != 0.0:
            label, grad, fitted, stangeby = ("beta_i", g_i, ALPHA_I,
                                             stangeby_beta(run))
        elif g_i == 0.0 and g_e != 0.0:
            label, grad, fitted, stangeby = ("alpha_e", g_e, ALPHA_E,
                                             ALPHA_E)
        elif math.isclose(g_i, g_e, rel_tol=1e-3):
            label, grad = "beta_i+alpha_e", g_i
            fitted = ALPHA_I + ALPHA_E
            stangeby = stangeby_beta(run) + ALPHA_E
        else:
            # Unequal gradients: one drift can't separate the two
            continue

        def coefficient(rate):
            return -run.m_imp * rate * run.drift / (run.charge**2 * grad)

        matches = [f for f in friction if same_plasma(f, run)]
        if matches:
            src = matches[0]
            nu, nu_err = src.rate, src.rate_err
            nu_from = src.name.replace("reg_", "")
        elif run.slow:
            nu, nu_err, nu_from = model_friction_rate(run), 0.0, "model"
        else:
            nu, nu_err, nu_from = fluid_friction_rate(run), 0.0, "fluid"
        flan = coefficient(nu)
        err = abs(flan) * math.hypot(run.drift_err / run.drift,
                                     nu_err / nu)
        # The model's expectation is for a slow impurity only (nan else)
        model = -run.m_imp * model_friction_rate(run) * run.drift_model \
            / (run.charge**2 * grad)
        rows.append(Coefficient(run, label,
                                coefficient(fluid_friction_rate(run)), flan,
                                err, nu_from, model, fitted, stangeby))
    return rows


def print_coefficients(rows):
    print("Effective thermal-force coefficients, F / (Z^2 grad T). "
          "Not checked: for comparison.")
    print("  Flan/fluid nu: with the fluid friction rate.  Flan/Flan nu: "
          "with Flan's own friction rate (nu from).")
    print("  model: what the collision model should give at the run's "
          "step (slow impurities only).")
    print("  fitted: heavy limit, what the skew is fitted to.  Stangeby: "
          "his mass-dependent beta_i (+ 0.71).")
    head = (f"{'run':<30} {'coefficient':<15} {'Flan/fluid nu':>13} "
            f"{'Flan/Flan nu':>12} {'+-':>6} {'nu from':<18} {'model':>6} "
            f"{'fitted':>6} {'Stangeby':>8}")
    print(head)
    print("-" * len(head))
    for r in rows:
        print(f"{r.run.name:<30} {r.label:<15} {r.flan_fluid:13.3f} "
              f"{r.flan:12.3f} {r.err:6.3f} {r.nu_from:<18} {r.model:6.3f} "
              f"{r.fitted:6.3f} {r.stangeby:8.3f}")


# ----------------------------------------------------------------------
# Plot
# ----------------------------------------------------------------------
COLORS = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300",
          "#7a5cd6", "#52514e"]
MARKERS = ["o", "s", "^", "D", "v", "P", "X", "h"]


def plot_friction(runs, path):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    ink, muted, grid, surface = "#0b0b0b", "#52514e", "#e4e3df", "#fcfcfb"
    plt.rcParams.update({"font.size": 10, "axes.edgecolor": muted,
                         "axes.labelcolor": ink, "xtick.color": muted,
                         "ytick.color": muted, "text.color": ink})
    fig, (ax, ax_r) = plt.subplots(
        2, 1, figsize=(8.5, 7.5), sharex=True, facecolor=surface,
        gridspec_kw={"height_ratios": [3, 1.4]})

    # Time in fluid slowing-down times and velocity in flow speeds, so
    # every variant should fall on the one fluid curve
    x_max = 0.0
    for k, run in enumerate(runs):
        color = COLORS[k % len(COLORS)]
        x = run.t[1:-1] * run.rate_fluid
        y = run.v_par[1:-1] / run.flow
        x_max = max(x_max, x[-1])
        step = max(1, len(x) // 40)
        label = (f"{run.name.replace('reg_friction_', '')}: W{run.charge}+, "
                 f"tau_s = {1e6 / run.rate:.2f} us "
                 f"(fluid {1e6 / run.rate_fluid:.2f})")
        ax.plot(x[::step], y[::step], MARKERS[k % len(MARKERS)], ms=5,
                color=color, mec=surface, mew=0.8, label=label, zorder=3)
        fluid = -np.expm1(-x)
        ax_r.plot(x[::step], (y - fluid)[::step], MARKERS[k % len(MARKERS)],
                  ms=4, color=color, mec=surface, mew=0.6, zorder=3)
        # What the collision model should give for this variant
        ax.plot(x, -np.expm1(-x * run.rate_model / run.rate_fluid),
                color=color, lw=1.0, alpha=0.8, zorder=2)

    x = np.linspace(0.0, x_max, 400)
    ax.plot(x, -np.expm1(-x), color=ink, lw=2.0, zorder=1,
            label="Fluid: 1 - exp(-t / tau_s)")
    ax_r.axhline(0.0, color=ink, lw=2.0, zorder=1)

    ax.set_ylabel("Mean parallel velocity / flow speed")
    ax.set_title("Friction: tungsten dragged up to the plasma flow\n"
                 "Markers: Flan.  Thin lines: expected from the collision "
                 "model at each run's time step.", loc="left", fontsize=11)
    ax.legend(frameon=False, loc="lower right")
    ax_r.set_ylabel("Flan - fluid")
    ax_r.set_xlabel("Time / fluid slowing-down time")
    for a in (ax, ax_r):
        a.set_facecolor(surface)
        a.grid(True, color=grid, lw=0.8)
        a.set_axisbelow(True)
        for side in ("top", "right"):
            a.spines[side].set_visible(False)
    fig.tight_layout()
    fig.savefig(path, dpi=150)
    plt.close(fig)


def plot_thermal_force(runs, path):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    ink, muted, grid, surface = "#0b0b0b", "#52514e", "#e4e3df", "#fcfcfb"
    plt.rcParams.update({"font.size": 10, "axes.edgecolor": muted,
                         "axes.labelcolor": ink, "xtick.color": muted,
                         "ytick.color": muted, "text.color": ink})
    fig, (ax, ax_b) = plt.subplots(
        1, 2, figsize=(13.0, 5.6), facecolor=surface,
        gridspec_kw={"width_ratios": [1.5, 1]})

    # Left: mean parallel velocity against time, with each variant's
    # fluid drift as a dashed line from where force balance is reached
    for k, run in enumerate(runs):
        color = COLORS[k % len(COLORS)]
        short = run.name.replace("reg_thermal_force_", "")
        x = run.t[1:-1] * run.rate_fluid
        step = max(1, len(x) // 60)
        ax.plot(x[::step], run.v_par[1:-1][::step],
                MARKERS[k % len(MARKERS)], ms=4, color=color, mec=surface,
                mew=0.6, label=short, zorder=3)
        ax.plot([5.0, x[-1]], [run.drift_fluid] * 2, color=color, lw=1.6,
                ls="--", zorder=2)
    ax.axvline(5.0, color=muted, lw=0.8, ls=":")
    ax.set_xlabel("Time / fluid slowing-down time")
    ax.set_ylabel("Mean parallel velocity (m/s)")
    ax.set_title("Thermal force: impurities pushed against the heat flux\n"
                 "Markers: Flan.  Dashed: fluid drift at force balance.",
                 loc="left", fontsize=11)
    ax.legend(frameon=False, ncol=2, loc="upper right")

    # Right: the drift at force balance, Flan against fluid and model
    y = np.arange(len(runs))[::-1]
    for k, run in enumerate(runs):
        color = COLORS[k % len(COLORS)]
        ax_b.errorbar(run.drift, y[k], xerr=run.drift_err,
                      fmt=MARKERS[k % len(MARKERS)], ms=7, color=color,
                      mec=surface, mew=0.8, capsize=3, zorder=3)
        ax_b.plot([run.drift_fluid] * 2, [y[k] - 0.32, y[k] + 0.32],
                  color=ink, lw=2.0, zorder=2)
        ax_b.plot([run.drift_model] * 2, [y[k] - 0.32, y[k] + 0.32],
                  color=color, lw=1.2, zorder=2)
    ax_b.set_yticks(y)
    ax_b.set_yticklabels([r.name.replace("reg_thermal_force_", "")
                          for r in runs])
    ax_b.set_xlabel("Drift at force balance (m/s)")
    ax_b.set_title("Markers: Flan.  Thick bar: fluid.\n"
                   "Thin bar: expected from the collision model.",
                   loc="left", fontsize=11)
    for a in (ax, ax_b):
        a.set_facecolor(surface)
        a.grid(True, color=grid, lw=0.8)
        a.set_axisbelow(True)
        for side in ("top", "right"):
            a.spines[side].set_visible(False)
    fig.tight_layout()
    fig.savefig(path, dpi=150)
    plt.close(fig)


# ----------------------------------------------------------------------
# Main
# ----------------------------------------------------------------------
def print_table(checks):
    head = (f"{'run':<30} {'check':<16} {'Flan':>10} {'expected':>10} "
            f"{'allowed':>9} {'unit':<4} result")
    print(head)
    print("-" * len(head))
    for c in checks:
        print(f"{c.run.name:<30} {c.what:<16} {c.value:10.3f} "
              f"{c.expected:10.3f} {c.tol:9.3f} {c.unit:<4} {c.status()}")


def main():
    directory = sys.argv[1] if len(sys.argv) > 1 else "."
    checks = []

    friction = [Run(p) for p in
                sorted(glob.glob(os.path.join(directory,
                                              "reg_friction_*.nc")))]
    for run in friction:
        checks += check_friction(run)
    plots = []
    if friction:
        plots.append(os.path.join(directory, "check_friction.png"))
        plot_friction(friction, plots[-1])

    thermal = [Run(p) for p in
               sorted(glob.glob(os.path.join(directory,
                                             "reg_thermal_force_*.nc")))]
    for run in thermal:
        checks += check_thermal_force(run)
    if thermal:
        plots.append(os.path.join(directory, "check_thermal_force.png"))
        plot_thermal_force(thermal, plots[-1])

    if not checks:
        print(f"No regression output files found in {directory}")
        return 1

    print_table(checks)
    print()
    coefficients = thermal_coefficients(thermal, friction)
    if coefficients:
        print_coefficients(coefficients)
        print()
    for path in plots:
        print(f"Plot: {path}")
    failed = [c for c in checks if c.enforced and not c.ok]
    missed = [c for c in checks if not c.enforced and not c.ok]
    print(f"\n{len(checks) - len(failed) - len(missed)} of {len(checks)} "
          f"checks passed, {len(missed)} expected misses, "
          f"{len(failed)} failed")
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main())
