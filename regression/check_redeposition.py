#!/usr/bin/env python3
"""
Checks the prompt redeposition regression cases (reg_prompt_redep_<variant>.nc,
made by reg_prompt_redep in regression/) against Fussmann's model and plots
them.

Sputtered tungsten atoms leave a wall with a cosine distribution, ionize
to W+ and gyrate in a field parallel to the wall, with no electric field
and no collisions. The fraction that comes back within its first gyration
depends only on p = e B / (m n_e S0), the mean ionization length over the
gyroradius, with S0 the ionization rate coefficient of W0:

    <r>(p) = (1 / p) int_0^inf exp(-x / p) / (1 + x^2) dx

(fussmann.r_exp, in this directory). Flan's prompt fraction is the rate of
prompt returns (within one W+ gyroperiod of first ionizing) over the
sputtered source.

Everything comes from the inputs saved in each output file: n_e, Te, B
and the impurity's mass from the slab options, and S0 from the same
OpenADAS file Flan read (openadas_root, openadas_year), interpolated the
same way (bilinearly in n_e and Te, clamped at the table's edges).

A run passes if its prompt fraction is within 4 standard errors (binomial,
from imp_num) plus 0.01 of <r>(p). The 0.01 allows for what the runs do
differently from the model on purpose: the field at 0.01 degrees to the
wall, and the time step. A run that lets atoms ionize past W+
(imp_iz_max_charge below 0, the "full" variant) is printed but not
checked, since the model leaves that out.

Runs with the sheath field on (sheath_model not "off"):

- drop0, with a drop of 1e-4 Te (sheath_drop), is Fussmann's case again and
  is checked the same way.
- <model>_a<angle> are printed, not checked: the prompt fraction, and the
  redeposited ions' mean charge, impact energy and angle, against the field's
  angle to the wall for each profile ("brooks", "borodkina"), plotted in
  check_redeposition_sheath.png. The difference between the profiles is the
  sensitivity to the profile's shape.
- brooks_ref, Brooks' reference case (Phys. Fluids B 2, 1858, 1990,
  Table 1), is printed beside his numbers, not checked. Flan's W0 ionization
  length there (n S0 from OpenADAS, and the sputtered atoms' speed) is about
  0.06 mm, several times shorter than Brooks' mean ionization height of
  0.40 mm, so its W ionizes deeper in the magnetic presheath: it gains less of
  the drop, comes back sooner, at a lower charge and nearer the normal. The
  difference is in the ionization rates, not the sheath.

Usage:
    python check_redeposition.py [directory with the .nc files]

Writes check_redeposition.png to that directory, prints a table, and exits
with status 1 if a check failed or no output files were found.
"""

import glob
import math
import os
import sys

import numpy as np
import netCDF4

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import fussmann  # noqa: E402

# The constants Flan uses (include/constants.h), so p is Flan's p
E_CHARGE = 1.602e-19   # C
AMU = 1.661e-27        # kg

N_SIGMA = 4.0   # statistical part of the tolerance
TOL_ABS = 0.01  # the rest

# ADF11 file names by atomic number (as OpenADAS::get_openadas_path)
SYMBOLS = {1: "h", 2: "he", 3: "li", 4: "be", 5: "b", 6: "c", 7: "n",
    8: "o", 9: "f", 10: "ne", 74: "w"}


# ----------------------------------------------------------------------
# OpenADAS
# ----------------------------------------------------------------------
def read_adf11(path):
    """
    Reads an OpenADAS ADF11 file (e.g., scd50_w.dat).

    Returns (ne, te, rates): ne (m^-3) and te (eV) of the table, and the
    rate coefficients (m^3/s) as an array (charge, te, ne), charge 0 being
    the first block (W0 -> W+ for scd).
    """
    with open(path) as f:
        lines = f.read().split("\n")
    head = lines[0].split()
    nd, nt, z1, z2 = int(head[1]), int(head[2]), int(head[3]), int(head[4])
    i = 2
    vals = []
    while len(vals) < nd + nt:
        vals += [float(x) for x in lines[i].split()]
        i += 1
    ne = 10.0 ** np.array(vals[:nd]) * 1e6
    te = 10.0 ** np.array(vals[nd:nd + nt])
    rates = []
    while len(rates) < z2 - z1 + 1:
        while "Z1=" not in lines[i]:
            i += 1
        i += 1
        block = []
        while len(block) < nd * nt:
            block += [float(x) for x in lines[i].split()]
            i += 1
        rates.append(10.0 ** np.array(block).reshape(nt, nd) * 1e-6)
    return ne, te, np.array(rates)


def rate_coeff(table, charge, ne, te):
    """
    A rate coefficient (m^3/s) at ne (m^-3) and te (eV), bilinear in ne and
    te and clamped at the table's edges, as OpenADAS::get_rate_coeff.
    """
    ne_t, te_t, rates = table

    def bracket(grid, v):
        v = min(max(v, grid[0]), grid[-1])
        k = int(np.searchsorted(grid, v, side="right")) - 1
        k = min(max(k, 0), len(grid) - 2)
        return k, (v - grid[k]) / (grid[k + 1] - grid[k])

    it, ft = bracket(te_t, te)
    idn, fn = bracket(ne_t, ne)
    r = rates[charge]
    return ((1 - ft) * (1 - fn) * r[it, idn] + ft * (1 - fn) * r[it + 1, idn]
        + (1 - ft) * fn * r[it, idn + 1] + ft * fn * r[it + 1, idn + 1])


# ----------------------------------------------------------------------
# One run
# ----------------------------------------------------------------------
class Run:
    """Inputs, p and the prompt fraction of one run."""

    def __init__(self, path):
        self.path = path
        self.name = os.path.splitext(os.path.basename(path))[0]
        with netCDF4.Dataset(path) as nc:
            inp = nc.groups["input"]

            def val(name):
                v = inp.variables[name][:]
                return v[0] if np.ndim(v) else v.item() if hasattr(v, "item") else v

            self.ne = float(val("slab_ne"))
            self.te = float(val("slab_te"))
            self.B = float(val("slab_bmag"))
            self.angle = float(val("slab_b_angle"))
            self.mass_amu = float(val("imp_mass_amu"))
            self.Z = int(val("imp_atom_num"))
            self.n = int(val("imp_num"))
            self.max_charge = int(val("imp_iz_max_charge"))
            names = inp.variables.keys()
            self.sheath = (str(val("sheath_model")) if "sheath_model" in names
                else "off")
            self.drop = (float(val("sheath_drop")) if "sheath_drop" in names
                else -1.0)
            self.adas_root = str(val("openadas_root"))
            self.adas_year = int(val("openadas_year"))
            sput = nc.groups["sputtering"]
            self.source = float(sput.getncattr("source_rate"))
            self.redep = float(sput.getncattr("redeposited_total"))
            self.prompt = float(sput.getncattr("prompt_total"))
            attrs = sput.ncattrs()
            self.mean_energy = (float(sput.getncattr("mean_impact_energy"))
                if "mean_impact_energy" in attrs else float("nan"))
            self.mean_angle = (float(sput.getncattr("mean_impact_angle"))
                if "mean_impact_angle" in attrs else float("nan"))
            self.mean_lifetime = (float(sput.getncattr("mean_lifetime"))
                if "mean_lifetime" in attrs else float("nan"))
            q = np.array(sput.variables["redeposited_charge"][:])
            self.mean_charge = (float((np.arange(len(q)) * q).sum() / q.sum())
                if q.sum() > 0 else float("nan"))

        year = "{:02d}".format(self.adas_year)
        adas = os.path.join(self.adas_root, "scd" + year,
            "scd{}_{}.dat".format(year, SYMBOLS[self.Z]))
        self.S0 = rate_coeff(read_adf11(adas), 0, self.ne, self.te)
        self.p = E_CHARGE * self.B / (self.mass_amu * AMU * self.ne * self.S0)
        self.expected = fussmann.r_exp(self.p)
        self.frac = self.prompt / self.source
        self.sigma = math.sqrt(max(self.frac * (1 - self.frac), 1e-12)
            / self.n)
        no_sheath = self.sheath == "off" or 0.0 < self.drop < 0.01
        self.checked = self.max_charge >= 0 and no_sheath
        self.scan = self.sheath != "off" and "_a" in self.name
        self.passed = abs(self.frac - self.expected) <= (N_SIGMA * self.sigma
            + TOL_ABS)


def main():
    directory = sys.argv[1] if len(sys.argv) > 1 else "."
    paths = sorted(glob.glob(os.path.join(directory,
        "reg_prompt_redep_*.nc")))
    if not paths:
        print("No reg_prompt_redep_*.nc files in " + directory)
        sys.exit(1)
    runs = sorted((Run(p) for p in paths), key=lambda r: r.p)

    print("{:<24} {:>9} {:>4} {:>9} {:>7} {:>9} {:>9} {:>8} {:>6}".format(
        "run", "ne (m-3)", "B", "S0 (m3/s)", "p", "prompt", "Fussmann",
        "tol", ""))
    failed = 0
    for r in runs:
        tol = N_SIGMA * r.sigma + TOL_ABS
        status = ("pass" if r.passed else "FAIL") if r.checked else "info"
        failed += r.checked and not r.passed
        print("{:<24} {:9.2e} {:4.1f} {:9.3e} {:7.3f} {:9.4f} {:9.4f} "
            "{:8.4f} {:>6}".format(r.name, r.ne, r.B, r.S0, r.p, r.frac,
            r.expected, tol, status))
    for r in runs:
        if r.redep > 0 and r.checked and r.prompt < 0.99 * r.redep:
            print("Note: in {}, {:.1%} of the returns weren't prompt".format(
                r.name, 1 - r.prompt / r.redep))
    for r in runs:
        if r.max_charge < 0 and r.sheath == "off":
            print("{}: atoms ionize past W+, {:.4f} back promptly against "
                "{:.4f} from the model".format(r.name, r.frac, r.expected))

    # The sheath angle scan
    scan = sorted((r for r in runs if r.scan), key=lambda r: (r.sheath,
        r.angle))
    if scan:
        print("\nWith the sheath field (not checked): Fussmann without it "
            "gives {:.4f}".format(scan[0].expected))
        print("{:<28} {:>6} {:>7} {:>7} {:>9} {:>9}".format("run",
            "B-wall", "prompt", "charge", "energy", "angle"))
        for r in scan:
            print("{:<28} {:6.1f} {:7.4f} {:7.2f} {:9.1f} {:9.1f}".format(
                r.name, r.angle, r.frac, r.mean_charge, r.mean_energy,
                r.mean_angle))

    # Brooks' reference case
    for r in runs:
        if r.name.endswith("brooks_ref"):
            print("\nBrooks' reference case (not checked; see the top of "
                "this script):")
            print("{:<26} {:>8} {:>8}".format("redeposited W", "Brooks",
                "Flan"))
            print("{:<26} {:8.1f} {:8.2f}".format("charge", 2.6,
                r.mean_charge))
            print("{:<26} {:8.0f} {:8.0f}".format("impact energy (eV)", 553,
                r.mean_energy))
            print("{:<26} {:8.0f} {:8.0f}".format("angle from normal (deg)",
                72, r.mean_angle))
            print("{:<26} {:8.2f} {:8.3f}".format("time to return (us)", 0.6,
                1e6 * r.mean_lifetime))

    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        fig, ax = plt.subplots(figsize=(5.5, 4))
        ps = np.logspace(-1.5, 1.5, 200)
        ax.plot(ps, [fussmann.r_exp(p) for p in ps], "k-",
            label=r"$\langle r\rangle(p)$, exponential average")
        ax.plot(ps, fussmann.r_fixed(ps), "k:",
            label=r"$1/(1+p^2)$, fixed length")
        chk = [r for r in runs if r.checked]
        ax.errorbar([r.p for r in chk], [r.frac for r in chk],
            yerr=[2 * r.sigma for r in chk], fmt="o", color="C0",
            label="Flan, ionizing to W+ only")
        info = [r for r in runs if r.max_charge < 0 and r.sheath == "off"]
        if info:
            ax.plot([r.p for r in info], [r.frac for r in info], "s",
                color="C3", label="Flan, full ionization")
        ax.set_xscale("log")
        ax.set_xlabel(r"$p = \lambda_{iz} / \rho$")
        ax.set_ylabel("fraction promptly redeposited")
        ax.set_ylim(0, 1.05)
        ax.legend(fontsize=8)
        fig.tight_layout()
        out = os.path.join(directory, "check_redeposition.png")
        fig.savefig(out, dpi=150)
        print("\nPlot: " + out)
    except ImportError:
        print("\nmatplotlib not found, no plot")

    if scan:
        try:
            import matplotlib
            matplotlib.use("Agg")
            import matplotlib.pyplot as plt
            fig, axs = plt.subplots(1, 3, figsize=(11, 3.5))
            for model, marker in (("brooks", "o-"), ("borodkina", "s--")):
                # The CPU runs; _gpu ones are in the table beside them
                rs = [r for r in scan if r.sheath == model
                      and not r.name.endswith("_gpu")]
                if not rs:
                    continue
                a = [r.angle for r in rs]
                axs[0].plot(a, [r.frac for r in rs], marker, label=model)
                axs[1].plot(a, [r.mean_energy for r in rs], marker,
                    label=model)
                axs[2].plot(a, [r.mean_angle for r in rs], marker,
                    label=model)
            axs[0].axhline(scan[0].expected, color="k", ls=":",
                label="Fussmann, no sheath")
            axs[0].set_ylabel("fraction promptly redeposited")
            axs[1].set_ylabel("mean impact energy (eV)")
            axs[2].set_ylabel("mean impact angle from normal (deg)")
            for ax in axs:
                ax.set_xlabel("angle between B and the wall (deg)")
            axs[0].legend(fontsize=8)
            fig.tight_layout()
            out = os.path.join(directory, "check_redeposition_sheath.png")
            fig.savefig(out, dpi=150)
            print("Plot: " + out)
        except ImportError:
            pass

    n_checked = sum(r.checked for r in runs)
    print("\n{} of {} checks passed, {} failed".format(n_checked - failed,
        n_checked, failed))
    sys.exit(1 if failed or n_checked == 0 else 0)


if __name__ == "__main__":
    main()
