#!/usr/bin/env python3
"""
Checks the sheath impact angle regression cases (reg_sheath_angles_<variant>.txt,
made by reg_sheath_angles in regression/) and plots them.

D+ test ions fall from 10 rho_s through Borodkina and Tsvetkov's sheath field
onto a wall, with the field at 45 or 60 degrees to the wall's normal. Two
checks per run:

1. Energy conservation (exact). Each ion should hit the wall with its starting
   energy plus e (phi(start) - phi_wall): B does no work and the field is
   static. The largest error over the ions must be below 1e-3 of the drop.
   This checks the field and the push near the wall, sub-stepping included.

2. The distribution of the angle from the normal, against Borodkina and
   Tsvetkov's model curves (Physics Procedia 71, 25, 2015, Fig. 2a, read off
   the figure by eye, to about 0.001 in dN/dtheta and 1 degree). The mean
   must be within 5 degrees and the standard deviation within 3 degrees.
   The paper gives the start only as "a shifted Maxwell velocity
   distribution" from 10 rho_s. Flan shifts it by c_s along B; shifting by
   0.5 c_s or 2 c_s instead moves the mean at 60 degrees from 31.5 to 29 or
   39 degrees. So the start alone could account for a few degrees, and the
   tolerance on the mean allows for it.

Usage:
    python check_sheath_angles.py [directory with the .txt files]

Writes check_sheath_angles.png to that directory, prints a table, and exits
with status 1 if a check failed or no output files were found.
"""

import glob
import os
import sys

import numpy as np

ENERGY_TOL = 1e-3   # of the drop
MEAN_TOL = 5.0      # degrees
STD_TOL = 3.0       # degrees

# Borodkina and Tsvetkov's model curves, (theta in degrees, dN/dtheta per
# degree), read off Fig. 2a
MODEL = {
    45.0: [(0, 0.002), (2, 0.005), (5, 0.008), (8, 0.0118), (12, 0.0145),
           (15, 0.017), (19, 0.0195), (23, 0.022), (27, 0.024), (31, 0.0255),
           (35, 0.0245), (39, 0.023), (42, 0.022), (45, 0.019), (48, 0.0165),
           (51, 0.015), (54, 0.012), (58, 0.0098), (61, 0.006), (65, 0.003),
           (68, 0.0015), (72, 0.0005), (76, 0.0), (90, 0.0)],
    60.0: [(0, 0.0), (5, 0.0015), (8, 0.003), (12, 0.0055), (15, 0.009),
           (18, 0.0125), (21, 0.016), (25, 0.0215), (28, 0.025), (31, 0.028),
           (34, 0.0295), (37, 0.0303), (40, 0.0295), (43, 0.027), (46, 0.024),
           (49, 0.021), (52, 0.016), (55, 0.011), (58, 0.0065), (61, 0.004),
           (64, 0.002), (68, 0.0008), (72, 0.0), (90, 0.0)],
}


def model_stats(alpha):
    """Mean and standard deviation (degrees) of a model curve, normalized."""
    t, f = np.array(MODEL[alpha]).T
    x = np.linspace(0.0, 90.0, 901)
    y = np.interp(x, t, f)
    n = np.trapz(y, x)
    mean = np.trapz(x * y, x) / n
    std = np.sqrt(np.trapz((x - mean) ** 2 * y, x) / n)
    return mean, std


class Run:
    """One reg_sheath_angles output file."""

    def __init__(self, path):
        self.name = os.path.splitext(os.path.basename(path))[0]
        self.info = {}
        bins = []
        with open(path) as f:
            for line in f:
                if line.startswith("#"):
                    continue
                key, val = line.split()
                try:
                    bins.append((int(key), float(val)))
                except ValueError:
                    self.info[key] = float(val)
        self.theta = np.array([b[0] for b in bins]) + 0.5
        counts = np.array([b[1] for b in bins])
        self.frac = counts / counts.sum()
        self.mean = float((self.theta * self.frac).sum())
        self.std = float(np.sqrt(((self.theta - self.mean) ** 2
            * self.frac).sum()))
        self.alpha = self.info["alpha"]
        self.model_mean, self.model_std = model_stats(self.alpha)
        self.energy_ok = self.info["max_energy_error"] <= ENERGY_TOL
        self.mean_ok = abs(self.mean - self.model_mean) <= MEAN_TOL
        self.std_ok = abs(self.std - self.model_std) <= STD_TOL


def main():
    directory = sys.argv[1] if len(sys.argv) > 1 else "."
    paths = sorted(glob.glob(os.path.join(directory,
        "reg_sheath_angles_*.txt")))
    if not paths:
        print("No reg_sheath_angles_*.txt files in " + directory)
        sys.exit(1)
    runs = [Run(p) for p in paths]

    print("{:<24} {:>5} {:>9} {:>11} {:>11} {:>6}".format("run", "alpha",
        "energy", "mean", "std", ""))
    print("{:<24} {:>5} {:>9} {:>11} {:>11}".format("", "(deg)",
        "error", "Flan/model", "Flan/model"))
    checks = 0
    failed = 0
    for r in runs:
        ok = r.energy_ok and r.mean_ok and r.std_ok
        checks += 3
        failed += (not r.energy_ok) + (not r.mean_ok) + (not r.std_ok)
        print("{:<24} {:5.0f} {:9.1e} {:5.1f}/{:<5.1f} {:5.1f}/{:<5.1f} "
            "{:>6}".format(r.name, r.alpha, r.info["max_energy_error"],
            r.mean, r.model_mean, r.std, r.model_std,
            "pass" if ok else "FAIL"))
        if r.info["hit"] < 0.99 * r.info["num_ions"]:
            print("Note: in {}, only {:.0f} of {:.0f} ions hit the wall".format(
                r.name, r.info["hit"], r.info["num_ions"]))

    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        fig, ax = plt.subplots(figsize=(5.5, 4))
        for r in runs:
            line, = ax.plot(r.theta, r.frac, "-",
                label="Flan, {:.0f}°".format(r.alpha))
            t, f = np.array(MODEL[r.alpha]).T
            ax.plot(t, f, "o--", ms=3, color=line.get_color(),
                label="Borodkina model, {:.0f}°".format(r.alpha))
        ax.set_xlabel("angle from the normal (degrees)")
        ax.set_ylabel(r"$dN/d\theta$ (per degree)")
        ax.set_xlim(0, 90)
        ax.legend(fontsize=8)
        fig.tight_layout()
        out = os.path.join(directory, "check_sheath_angles.png")
        fig.savefig(out, dpi=150)
        print("\nPlot: " + out)
    except ImportError:
        print("\nmatplotlib not found, no plot")

    print("\n{} of {} checks passed, {} failed".format(checks - failed, checks,
        failed))
    sys.exit(1 if failed else 0)


if __name__ == "__main__":
    main()
