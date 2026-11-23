#!/usr/bin/env python3
"""
Checks the self-sputtering regression cases (reg_self_sputter_<variant>.nc,
made by reg_self_sputter in regression/) and plots them.

Every W atom or ion that lands on the slab wall sputters Y W atoms, and
the particle carries on as one of them with its weight times Y. The output
has the weight carried into each generation g, W_g (W_0 = S, the D+
sputtered source), and the weight landing on the wall from it, L_g. The
return fraction of generation g is f_g = L_g / W_g.

With a constant yield Y0 (y05, y125), every particle of a generation has
the same weight: w_0 = S / imp_num, then w_{g+1} = w_g Y0, or w_min (the
roulette weight, sput_self_min_weight w_0) where w_g Y0 is below w_min.
So W_g / w_g counts the particles in generation g, and the checks are:

1. Bookkeeping: W_{g+1} = Y0 L_g, exactly where no roulette is played,
   and otherwise the count in g+1 is binomial, from the n_g = L_g / w_g
   landings each carrying on with probability w_g Y0 / w_min (4 standard
   deviations plus 1). The counts have to come out whole numbers.
2. The tallies agree: the self-sputtered total and its sum over the
   elements are Y0 sum_g L_g, the redeposited total is sum_g L_g.
3. The self-sputtered atoms come back as often as the D-sputtered ones:
   f_0 against f_1, the pooled f_g for g >= 1, within 4 standard errors
   (binomial, from the counts) plus TOL_F. Prompt redeposition doesn't
   depend on the atom's speed (p = Omega tau_iz), but the W-sputtered
   atoms are faster, so slightly more of them are still out in the plasma
   when they're stopped (imp_max_time). TOL_F allows for that.
4. The chain adds up: the self-sputtered rate T = sum_{g>=1} W_g against
   S f_0 Y0 / (1 - f_1 Y0), and the multiplication T / (S + T) against
   f_1 Y0, to TOL_T.

The RustBCA table runs (w_te<Te>) are printed, not checked: the
self-sputtered rate and multiplication, the mean yield and impact energy,
and the yield at the energies the landing charges fall through,
sum_q P(q) Y(q phi_fl + 2 Ti, mean angle), with phi_fl the floating drop
for the slab's D plasma.

Usage:
    python check_self_sputter.py [directory with the .nc files]

Writes check_self_sputter.png to that directory, prints a table, and exits
with status 1 if a check failed or no constant-yield outputs were found.
"""

import glob
import math
import os
import sys

import h5py
import numpy as np
import netCDF4

N_SIGMA = 4.0   # statistical part of the tolerances
TOL_F = 0.02    # f_0 against f_1, beyond the statistics
TOL_T = 0.02    # relative, for the chain's total and multiplication
TOL_EXACT = 1e-9  # relative, for the sums that should agree exactly

ME_AMU = 0.000548


def interp_yield(path, energy, angle):
    """
    Interpolates a yield table as Flan does (YieldTable::get_yield):
    bilinear in log10(E) and angle, 0 below the lowest energy, the highest
    energy above the table, the angle clamped to the table's.
    """
    with h5py.File(path, "r") as f:
        e = np.array(f["energy"])
        a = np.array(f["angle"])
        y = np.array(f["yield"])
    if not energy >= e[0]:
        return 0.0
    le = min(math.log10(energy), math.log10(e[-1]))
    lgrid = np.log10(e)
    angle = min(max(angle, a[0]), a[-1])
    i = min(max(np.searchsorted(lgrid, le, side="right") - 1, 0),
        len(e) - 2)
    j = min(max(np.searchsorted(a, angle, side="right") - 1, 0),
        len(a) - 2)
    fe = (le - lgrid[i]) / (lgrid[i + 1] - lgrid[i])
    fa = (angle - a[j]) / (a[j + 1] - a[j])
    return ((1 - fe) * ((1 - fa) * y[i, j] + fa * y[i, j + 1])
        + fe * ((1 - fa) * y[i + 1, j] + fa * y[i + 1, j + 1]))


class Run:
    """One reg_self_sputter output and its checks."""

    def __init__(self, path):
        self.path = path
        self.name = os.path.splitext(os.path.basename(path))[0] \
            .replace("reg_self_sputter_", "")
        self.directory = os.path.dirname(path)
        with netCDF4.Dataset(path) as nc:
            inp = nc.groups["input"]

            def val(name):
                v = inp.variables[name][:]
                return v[0] if np.ndim(v) else v.item() if hasattr(v, "item") \
                    else v

            self.n = int(val("imp_num"))
            self.te = float(val("slab_te"))
            self.ti = float(val("slab_ti"))
            self.mi_amu = float(val("gkyl_ion_mass_amu"))
            self.table = str(val("sput_self_table_file"))
            sput = nc.groups["sputtering"]
            att = sput.getncattr
            self.S = float(att("source_rate"))
            self.T = float(att("self_sputtered_total"))
            self.w_min = float(att("self_min_weight"))
            self.landed = float(att("self_landed"))
            self.redep = float(att("redeposited_total"))
            self.mean_yield = float(att("self_mean_yield"))
            self.mean_energy = float(att("self_mean_impact_energy"))
            self.mean_angle = float(att("self_mean_impact_angle"))
            self.mult = float(att("self_multiplication"))
            self.W = np.array(sput.variables["generation_weight"][:])
            self.L = np.array(sput.variables["generation_landed"][:])
            self.per_elem = np.array(sput.variables["self_sputtered"][:])
            q = np.array(sput.variables["redeposited_charge"][:])
            self.charge = q / q.sum() if q.sum() > 0 else q

        # A constant-yield table?
        tpath = self.table if os.path.isabs(self.table) \
            else os.path.join(self.directory, self.table)
        self.table_path = tpath
        with h5py.File(tpath, "r") as f:
            y = np.array(f["yield"])
        self.constant = bool(np.all(y == y.flat[0]))
        self.Y0 = float(y.flat[0])
        self.checks = []
        if self.constant:
            self.check_constant()

    def add(self, name, ok, detail):
        self.checks.append((name, bool(ok), detail))

    def check_constant(self):
        """Checks 1 to 4 of the module docstring."""
        Y0, W, L = self.Y0, self.W.copy(), self.L
        W[0] = self.S
        G = int(np.max(np.nonzero(W)[0]))

        # Each generation's particle weight
        w = np.zeros(len(W))
        w[0] = self.S / self.n
        for g in range(1, len(W)):
            w[g] = w[g - 1] * Y0
            if w[g] < self.w_min * (1 - 1e-12):
                w[g] = self.w_min
        N = W / w
        nl = L / w
        self.N, self.nl = N, nl

        # 1. Whole counts, and W_{g+1} from L_g
        frac = max(np.max(np.abs(N[:G + 1] - np.round(N[:G + 1]))),
            np.max(np.abs(nl[:G + 1] - np.round(nl[:G + 1]))))
        self.add("whole counts", frac < 1e-6, "largest fraction {:.1e}"
            .format(frac))
        worst = 0.0
        ok = True
        for g in range(G + 1):
            if g + 1 >= len(W):
                break
            p = w[g] * Y0 / w[g + 1]
            if p >= 1 - 1e-12:
                err = abs(W[g + 1] - Y0 * L[g]) / max(Y0 * L[g], 1e-300)
                ok &= err < TOL_EXACT if L[g] > 0 else W[g + 1] == 0
                worst = max(worst, err / TOL_EXACT)
            else:
                n_land = round(nl[g])
                mean = n_land * p
                sd = math.sqrt(n_land * p * (1 - p))
                dev = abs(round(N[g + 1]) - mean)
                ok &= dev <= N_SIGMA * sd + 1
                worst = max(worst, dev / (N_SIGMA * sd + 1))
        self.add("W_{g+1} = Y0 L_g", ok, "worst {:.2f} of the tolerance"
            .format(worst))

        # 2. The tallies agree
        sumL = L.sum()
        errs = [abs(self.T - Y0 * sumL) / (Y0 * sumL),
            abs(self.per_elem.sum() - self.T) / self.T,
            abs(self.redep - sumL) / sumL,
            abs(self.landed - sumL) / sumL]
        self.add("tallies agree", max(errs) < TOL_EXACT, "largest {:.1e}"
            .format(max(errs)))

        # 3. f_0 against f_1
        self.f0 = L[0] / W[0]
        self.f1 = L[1:].sum() / W[1:].sum()
        n0, n1 = N[0], N[1:].sum()
        f = 0.5 * (self.f0 + self.f1)
        sd = math.sqrt(f * (1 - f) * (1 / n0 + 1 / n1))
        self.add("f_1 = f_0", abs(self.f1 - self.f0) <= N_SIGMA * sd + TOL_F,
            "f_0 {:.4f}, f_1 {:.4f} (+- {:.4f})".format(self.f0, self.f1,
            sd))

        # 4. The chain's total and multiplication
        self.T_pred = self.S * self.f0 * Y0 / (1 - self.f1 * Y0)
        self.T_model = self.S * self.f0 * Y0 / (1 - self.f0 * Y0)
        err = abs(self.T - self.T_pred) / self.T_pred
        self.add("T = S f0 Y0 / (1 - f1 Y0)", err < TOL_T,
            "T/S {:.4f}, predicted {:.4f} ({:.4f} with f_0 alone)".format(
            self.T / self.S, self.T_pred / self.S, self.T_model / self.S))
        err_m = abs(self.mult - self.f1 * Y0) / (self.f1 * Y0)
        self.add("M = f1 Y0", err_m < TOL_T, "M {:.4f}, f1 Y0 {:.4f}".format(
            self.mult, self.f1 * Y0))

    def estimated_yield(self):
        """sum_q P(q) Y(q phi_fl + 2 Ti, mean angle), for the table runs."""
        lam = -0.5 * math.log(2 * math.pi * ME_AMU / self.mi_amu
            * (1 + self.ti / self.te))
        phi = lam * self.te
        y = 0.0
        for q, p in enumerate(self.charge):
            if p > 0:
                y += p * interp_yield(self.table_path, q * phi + 2 * self.ti,
                    self.mean_angle)
        return y, phi

    @property
    def passed(self):
        return all(ok for _, ok, _ in self.checks)


def plot(runs, directory):
    """Weight per generation over the source, with (f Y0)^g for the
    constant-yield runs."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(figsize=(6, 4))
    for r in runs:
        W = r.W.copy()
        W[0] = r.S
        g = np.nonzero(W)[0]
        line, = ax.semilogy(g, W[g] / r.S, "o", label=r.name)
        if r.constant:
            gg = np.arange(0, g.max() + 1)
            model = np.where(gg == 0, 1.0, r.f0 * r.Y0
                * (r.f1 * r.Y0) ** np.maximum(gg - 1, 0))
            ax.semilogy(gg, model, "-", color=line.get_color(), lw=1)
    ax.set_xlabel("generation g")
    ax.set_ylabel("W_g / S")
    ax.set_title("Self-sputtering: weight per generation\n"
        "(lines: f0 Y0 (f1 Y0)^(g-1))", fontsize=10)
    ax.legend(fontsize=8)
    fig.tight_layout()
    out = os.path.join(directory, "check_self_sputter.png")
    fig.savefig(out, dpi=120)
    return out


def main():
    directory = sys.argv[1] if len(sys.argv) > 1 else "."
    paths = sorted(glob.glob(os.path.join(directory,
        "reg_self_sputter_*.nc")))
    if not paths:
        print("No reg_self_sputter_*.nc files in " + directory)
        sys.exit(1)
    runs = [Run(p) for p in paths]

    failed = False
    n_checked = 0
    for r in runs:
        if r.constant:
            n_checked += 1
            print("{}: Y0 = {}, S = {:.4g} atoms/s, T/S = {:.4f}, "
                "M = {:.4f}, deepest generation {}".format(r.name, r.Y0, r.S,
                r.T / r.S, r.mult, int(np.max(np.nonzero(r.W)[0]))))
            for name, ok, detail in r.checks:
                print("  {:4s} {:28s} {}".format("pass" if ok else "FAIL",
                    name, detail))
                failed |= not ok
        else:
            y_est, phi = r.estimated_yield()
            mean_q = float((np.arange(len(r.charge)) * r.charge).sum())
            print("{} (info): Te = {} eV, S = {:.4g} atoms/s, T/S = {:.4f}, "
                "M = {:.4f}".format(r.name, r.te, r.S, r.T / r.S, r.mult))
            print("  landing: mean charge {:.2f}, impact energy {:.1f} eV, "
                "angle {:.1f} deg, mean yield {:.4f}".format(mean_q,
                r.mean_energy, r.mean_angle, r.mean_yield))
            print("  sum_q P(q) Y(q phi_fl + 2 Ti) = {:.4f} (phi_fl = {:.1f} "
                "V)".format(y_est, phi))

    out = plot(runs, directory)
    print("Wrote " + out)
    if n_checked == 0:
        print("No constant-yield runs (y05, y125) to check.")
    print("FAILED" if failed or n_checked == 0 else "All checks passed")
    sys.exit(1 if failed or n_checked == 0 else 0)


if __name__ == "__main__":
    main()
