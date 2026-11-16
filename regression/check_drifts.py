#!/usr/bin/env python3
"""
Checks the drift regression cases (reg_drifts_<variant>.nc, made by
reg_drifts in regression/) against single-particle theory and plots them.

Every particle of a run starts at the same place with the same velocity
and there are no collisions, so the mean velocity in the output is the
velocity of one particle: a steady drift, plus the gyration around it.
The two are separated by a least-squares fit of each velocity component
to

    v(t) = c + d t + gyration at the cyclotron frequency and twice it,

with the frequency itself fitted near q B / m. c is the drift (at t = 0)
and d its rate of change.

What is expected follows from the inputs saved in each output file, not
from the name of the variant. For a positive ion of charge q and mass m:

  Straight slab, B along +X
    gyration       turns +Y toward -Z at omega = q B / m, speed unchanged
    E along +Y     ExB drift of E / B along -Z
    dE/dt          polarization drift of (m / (q B^2)) dE/dt along +Y
    dB/dY          grad-B drift of m v_perp^2 g / (2 q B) along +Z, with
                   g = (dB/dY) / B
  Bent slab, B along the toroidal angle and falling as 1 / R
    curvature      drift of m v_par^2 / (q B R) along +Z (upward), plus
                   the grad-B drift m v_perp^2 / (2 q B R) the same way

Usage:
    python check_drifts.py [directory with the .nc files]

Writes check_drifts.png to that directory, prints a table, and exits with
status 1 if a check failed or no output files were found.
"""

import glob
import math
import os
import sys

import numpy as np
import netCDF4

E_CHARGE = 1.602176634e-19   # C
AMU = 1.66053906660e-27      # kg

TOL_DRIFT = 0.02   # drifts, relative to the size of the run's main drift
TOL_GYRO = 0.005   # cyclotron frequency and speed in the gyration case


# ----------------------------------------------------------------------
# Reading a Flan output file
# ----------------------------------------------------------------------
class Run:
    """Inputs and the mean impurity velocity against time of one run."""

    def __init__(self, path):
        self.path = path
        self.name = os.path.splitext(os.path.basename(path))[0]
        self.short = self.name.replace("reg_drifts_", "")
        with netCDF4.Dataset(path) as nc:
            inp = nc["input"]

            def scalar(key):
                return float(np.ravel(inp[key][:])[0])

            self.bmag = scalar("slab_bmag")                # T
            self.e0 = scalar("slab_efield_y")              # V/m
            self.e_rate = scalar("slab_efield_y_rate")     # V/m/s
            self.bgrad = scalar("slab_bgrad_y")            # 1/m
            self.r0 = scalar("slab_major_radius")          # m, 0 = straight
            self.mass = scalar("imp_mass_amu") * AMU       # kg
            self.q = scalar("imp_init_charge") * E_CHARGE  # C
            self.v0 = np.array([scalar("imp_vX_start_val"),
                                scalar("imp_vY_start_val"),
                                scalar("imp_vZ_start_val")])   # m/s
            self.start = np.array([scalar("imp_xstart_val"),
                                   scalar("imp_ystart_val"),
                                   scalar("imp_zstart_val")])

            # Mean velocity of all the particles in each time frame: the
            # cell averages weighted by the counts in each cell. v_X, v_Y,
            # v_Z are Cartesian. v_x, v_y, v_z are the rates of change of
            # the computational coordinates, which for the bent slab are
            # R (m/s), Z (m/s) and the toroidal angle (rad/s).
            self.t = np.array(nc["geometry"]["time"][:], dtype=float)
            counts = np.array(nc["output"]["Nz"][:], dtype=float)
            total = counts.sum(axis=(1, 2, 3))
            if np.any(total[1:-1] <= 0):
                raise RuntimeError(f"{self.name}: particles were lost "
                                   "before the end of the run")
            total = np.where(total > 0, total, np.nan)
            self.v = {}
            for key in ("v_X", "v_Y", "v_Z", "v_x", "v_y", "v_z"):
                v = np.array(nc["output"][key][:], dtype=float)
                self.v[key] = (v * counts).sum(axis=(1, 2, 3)) / total

    @property
    def bent(self):
        return self.r0 > 0.0

    def b_start(self):
        """Field strength where the particles start (T)."""
        if self.bent:
            return self.bmag * self.r0 / self.start[0]   # x is R
        return self.bmag * (1.0 + self.bgrad * self.start[1])

    def omega(self):
        """Cyclotron frequency q B / m where the particles start."""
        return self.q * self.b_start() / self.mass


# ----------------------------------------------------------------------
# Separating drift from gyration
# ----------------------------------------------------------------------
def fit(run, keys):
    """
    Fits the velocity components `keys` to c + d t + harmonics of the
    gyration (see the top of the file), over the frames that cover a full
    frame width (all but the first and last). The gyration frequency is
    the one that fits all the components best, searched within 5% of
    q B / m. Returns {key: (c, d, amplitude of the gyration)} and the
    frequency.

    Each frame is an average over its width, which doesn't change c or d
    and lowers the gyration's amplitude by sin(x) / x, x = omega * width
    / 2. The amplitude returned is corrected for that.
    """
    t = run.t[1:-1]

    def solve(omega):
        basis = np.column_stack([
            np.ones_like(t), t,
            np.cos(omega * t), np.sin(omega * t),
            np.cos(2.0 * omega * t), np.sin(2.0 * omega * t)])
        out, resid = {}, 0.0
        for key in keys:
            coef, *_ = np.linalg.lstsq(basis, run.v[key][1:-1], rcond=None)
            resid += np.sum((basis @ coef - run.v[key][1:-1])**2)
            out[key] = coef
        return resid, out

    omega0 = run.omega()
    grid = omega0 * np.linspace(0.95, 1.05, 2001)
    resids = [solve(w)[0] for w in grid]
    omega = grid[int(np.argmin(resids))]
    _, coefs = solve(omega)

    x = 0.5 * omega * (run.t[1] - run.t[0])
    smear = math.sin(x) / x
    result = {key: (c[0], c[1], math.hypot(c[2], c[3]) / smear, c[2] / smear,
                    c[3] / smear) for key, c in coefs.items()}
    return result, omega


# ----------------------------------------------------------------------
# Checks
# ----------------------------------------------------------------------
class Check:
    """One line of the results table."""

    def __init__(self, run, what, value, expected, tol, unit):
        self.run, self.what, self.unit = run, what, unit
        self.value, self.expected, self.tol = value, expected, tol
        self.ok = abs(value - expected) <= tol


def check_straight(run):
    """Straight slab: B along +X, so the motion of interest is in Y, Z."""
    fits, omega = fit(run, ("v_Y", "v_Z"))
    run.fits, run.omega_fit = fits, omega
    b, q, m = run.b_start(), run.q, run.mass
    checks = []

    # Expected drifts. The speed of gyration that enters the grad-B drift
    # is the starting speed across the field (there is no E field in the
    # grad-B case, so no drift to take off first).
    v_perp_sq = run.v0[1]**2 + run.v0[2]**2
    exb = -run.e0 / b                                 # along Z, at t = 0
    exb_rate = -run.e_rate / b                        # its rate of change
    polarization = m * run.e_rate / (q * b * b)       # along Y
    gradb = m * v_perp_sq * run.bgrad / (2.0 * q * b)  # along Z
    run.expected = {"v_Y": (polarization, 0.0),
                    "v_Z": (exb + gradb, exb_rate)}

    drifting = max(abs(exb), abs(polarization), abs(gradb)) > 0.0
    if not drifting:
        # Pure gyration: frequency, speed, and the sense of rotation.
        # Starting along +Y, v_Y = v cos(omega t) and v_Z = -v sin(omega t)
        speed = math.sqrt(v_perp_sq)
        checks.append(Check(run, "omega", omega * 1e-6, run.omega() * 1e-6,
                            TOL_GYRO * run.omega() * 1e-6, "rad/us"))
        checks.append(Check(run, "speed (v_Y)", fits["v_Y"][2], speed,
                            TOL_GYRO * speed, "m/s"))
        checks.append(Check(run, "speed (v_Z)", fits["v_Z"][2], speed,
                            TOL_GYRO * speed, "m/s"))
        along = (run.v0[1] * fits["v_Y"][3] + run.v0[2] * fits["v_Z"][3]) \
            / speed
        turn = (run.v0[1] * fits["v_Z"][4] - run.v0[2] * fits["v_Y"][4]) \
            / speed
        checks.append(Check(run, "v along start", along, speed,
                            TOL_GYRO * speed, "m/s"))
        checks.append(Check(run, "turn, +Y to -Z", turn, -speed,
                            TOL_GYRO * speed, "m/s"))
        scale = speed * TOL_GYRO / TOL_DRIFT
    else:
        scale = max(abs(exb), abs(polarization), abs(gradb))

    # Drifts, and no motion along the field
    tol = TOL_DRIFT * scale
    checks.append(Check(run, "drift along Y", fits["v_Y"][0], polarization,
                        tol, "m/s"))
    checks.append(Check(run, "drift along Z", fits["v_Z"][0], exb + gradb,
                        tol, "m/s"))
    if run.e_rate != 0.0:
        # Rate of change of the ExB drift, in m/s per microsecond
        checks.append(Check(run, "d(drift Z)/dt", fits["v_Z"][1] * 1e-6,
                            exb_rate * 1e-6,
                            TOL_DRIFT * abs(exb_rate) * 1e-6, "m/s/us"))
    checks.append(Check(run, "mean v_X", np.mean(run.v["v_X"][1:-1]),
                        run.v0[0], tol, "m/s"))
    return checks


def check_bent(run):
    """
    Bent slab: B along the toroidal angle. The drift is in Z, which is
    the computational y, and v_x is the radial velocity.
    """
    fits, omega = fit(run, ("v_x", "v_y"))
    run.fits, run.omega_fit = fits, omega
    b, q, m = run.b_start(), run.q, run.mass
    r, phi = run.start[0], run.start[2]

    # Starting velocity along the field (the toroidal direction, (-sin
    # phi, cos phi, 0)) and across it
    v_par = -run.v0[0] * math.sin(phi) + run.v0[1] * math.cos(phi)
    v_perp_sq = float(run.v0 @ run.v0) - v_par**2
    drift = m * (v_par**2 + 0.5 * v_perp_sq) / (q * b * r)
    run.expected = {"v_x": (0.0, 0.0), "v_y": (drift, 0.0)}

    tol = TOL_DRIFT * abs(drift)
    checks = [
        Check(run, "drift along Z", fits["v_y"][0], drift, tol, "m/s"),
        Check(run, "drift along R", fits["v_x"][0], 0.0, tol, "m/s"),
        # The particle keeps going around at its parallel speed: v_z is
        # the rate of change of the toroidal angle
        Check(run, "v_par (R dphi/dt)", r * np.mean(run.v["v_z"][1:-1]),
              v_par, TOL_GYRO * abs(v_par), "m/s"),
    ]
    return checks


# ----------------------------------------------------------------------
# Plot
# ----------------------------------------------------------------------
def plot(runs, path):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    ink, muted, grid, surface = "#0b0b0b", "#52514e", "#e4e3df", "#fcfcfb"
    colors = ["#2a78d6", "#eb6834"]
    plt.rcParams.update({"font.size": 10, "axes.edgecolor": muted,
                         "axes.labelcolor": ink, "xtick.color": muted,
                         "ytick.color": muted, "text.color": ink})
    ncol = min(3, len(runs))
    nrow = math.ceil(len(runs) / ncol)
    fig, axes = plt.subplots(nrow, ncol, figsize=(5.2 * ncol, 4.6 * nrow),
                             facecolor=surface, squeeze=False)
    labels = {"v_Y": "v_Y", "v_Z": "v_Z", "v_x": "v_R", "v_y": "v_Z"}

    for ax, run in zip(axes.flat, runs):
        t_us = run.t[1:-1] * 1e6
        for k, key in enumerate(run.fits):
            c, d = run.fits[key][0], run.fits[key][1]
            c_exp, d_exp = run.expected[key]
            ax.plot(t_us, run.v[key][1:-1], "o", ms=2.5, color=colors[k],
                    alpha=0.55, zorder=2)
            ax.plot(t_us, c + d * run.t[1:-1], color=colors[k], lw=2.0,
                    zorder=3, label=f"{labels[key]}: drift {c:.1f} m/s "
                    f"(theory {c_exp:.1f})")
            ax.plot(t_us, c_exp + d_exp * run.t[1:-1], color=ink, lw=1.0,
                    ls="--", zorder=4)
        ax.set_title(run.short, loc="left", fontsize=11)
        ax.set_xlabel("Time (us)")
        ax.set_ylabel("Mean velocity (m/s)")
        ax.legend(frameon=False, fontsize=8, loc="upper center",
                  bbox_to_anchor=(0.5, -0.17), ncol=1)
        ax.set_facecolor(surface)
        ax.grid(True, color=grid, lw=0.8)
        ax.set_axisbelow(True)
        for side in ("top", "right"):
            ax.spines[side].set_visible(False)
    for ax in axes.flat[len(runs):]:
        ax.set_visible(False)
    fig.suptitle("Drifts. Dots: Flan (drift plus gyration).  Solid: fitted "
                 "drift.  Dashed: theory.", x=0.02, ha="left", fontsize=12)
    fig.tight_layout(rect=(0.0, 0.0, 1.0, 0.95))
    fig.savefig(path, dpi=150)
    plt.close(fig)


# ----------------------------------------------------------------------
# Main
# ----------------------------------------------------------------------
def main():
    directory = sys.argv[1] if len(sys.argv) > 1 else "."
    runs = [Run(p) for p in
            sorted(glob.glob(os.path.join(directory, "reg_drifts_*.nc")))]
    if not runs:
        print(f"No reg_drifts output files found in {directory}")
        return 1

    checks = []
    for run in runs:
        checks += check_bent(run) if run.bent else check_straight(run)

    head = (f"{'run':<24} {'check':<18} {'Flan':>11} {'expected':>11} "
            f"{'allowed':>9} {'unit':<6} result")
    print(head)
    print("-" * len(head))
    for c in checks:
        print(f"{c.run.name:<24} {c.what:<18} {c.value:11.3f} "
              f"{c.expected:11.3f} {c.tol:9.3f} {c.unit:<6} "
              f"{'pass' if c.ok else 'FAIL'}")

    path = os.path.join(directory, "check_drifts.png")
    plot(runs, path)
    print(f"\nPlot: {path}")
    failed = [c for c in checks if not c.ok]
    print(f"\n{len(checks) - len(failed)} of {len(checks)} checks passed, "
          f"{len(failed)} failed")
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main())
