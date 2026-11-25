"""
Checks the reg_sheath_impacts runs: the sputtered source from the D+
impact energies and angles of the sheath model (sput_impact_opt =
"sheath"). Usage:

	python check_sheath_impacts.py <directory>

It reads reg_sheath_impacts_*.nc in the directory (and, for the bor
variants, reg_sheath_angles_a45_sampler.txt and _a60_sampler.txt from
reg_sheath_angles), takes the inputs from the files, prints a table and
writes check_sheath_impacts.png there. The exit status is 1 if a check
failed.

What the test ions start with (Sheath::sample_impacts): 10 rho_s from the
wall, a Maxwellian at Ti shifted by c_s = sqrt(e (Te + Ti) / m_i) along B,
keeping those moving toward the wall along B, each weighted by its speed
along B (its share of the flux). The field is static and B does no work,
so each ion hits with its starting energy plus e (phi(start) - phi_wall).
With x the speed along B over v_th = sqrt(e Ti / m_i) and s = c_s / v_th,
the flux-weighted mean starting energy is

	E_start = Ti + (Ti / 2) M3 / M1,   M_k = int_0^inf x^k exp(-(x - s)^2 / 2) dx,

the two directions across B giving Ti, and its spread follows the same
way. The checks:

 - const (Y = 0.02): "sheath" and "formula" give the same source, and the
   yield is 0.02.
 - lin (Y = 1e-3 E / eV, flat in angle): the yield is 1e-3 times the mean
   impact energy, and that is the drop plus E_start, within 4 standard
   errors (from the spread of the starting energies and the number of
   test ions) plus 1e-3 of the drop for the push.
 - lin_a0: with B along the normal the field only speeds the ions up
   along it, so tan(angle) = v_perp / sqrt(v_par^2 + 2 e drop / m_i); the
   flux-weighted mean angle is worked out from that by Monte Carlo.
 - bor: the mean angle (weighted by the flux) matches the
   one-dimensional sampler at the same plasma (reg_sheath_angles), within
   4 standard errors and GRID_TOL for the impact table's grid in
   lambda_D / rho_s (a factor of 2), which these cases are sensitive to
   since their Debye length is large.
 - Every "sheath" run: no ion lost, energy conserved to 1e-3 of the drop,
   and the histograms add up to the ion flux and the source.

The dw runs (RustBCA's D on W) are printed, "formula" beside "sheath",
not checked.
"""

import argparse
import glob
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import netCDF4
import numpy as np


LIN_SLOPE = 1e-3        # make_impact_tables.py
CONST_YIELD = 0.02
START_RHO_S = 10.0      # Sheath::sample_impacts' default
ENERGY_TOL = 1e-3       # of the drop
SAME_TOL = 1e-9         # relative, for the same source
INTERP_TOL = 1e-4       # relative, linear yield interpolated in log E
HIST_TOL = 1e-9         # relative, histogram sums
# Degrees, on top of 4 standard errors, for bor: the table takes
# lambda_D / rho_s to the nearest power of 2, up to a factor sqrt(2) off.
# In the 45 degree case (lambda_D / rho_s = 0.18) a factor sqrt(2) moves
# the mean angle by 1.5 to 1.9 degrees, in the 60 degree case (0.011) by
# 0.3 (the sampler with 100000 ions).
GRID_TOL = 2.0
NUM_SE = 4.0


class Run:
	"""
	One reg_sheath_impacts output file: the plasma, the options and what
	the sputtered source came out as.
	"""

	def __init__(self, path):
		self.path = path
		self.name = os.path.basename(path)[len("reg_sheath_impacts_"):-3]
		with netCDF4.Dataset(path) as nc:
			inp = nc["input"]

			def num(key):
				return float(np.ravel(inp[key][...])[0])

			def text(key):
				return str(np.ravel(inp[key][...])[0])

			self.te = num("slab_te")
			self.ti = num("slab_ti")
			self.theta = 90.0 - num("slab_b_angle")
			self.mi = num("gkyl_ion_mass_amu")
			self.me = num("gkyl_elec_mass_amu")
			self.model = text("sheath_model")
			self.opt = text("sput_impact_opt")

			g = nc["sputtering"]
			self.attrs = {k: g.getncattr(k) for k in g.ncattrs()}
			self.source = float(self.attrs["source_rate"])
			flux = np.asarray(g["flux"][:], dtype=float)
			area = np.asarray(g["area"][:], dtype=float)
			w = flux * area[None, :]
			self.ions = w.sum(axis=1).mean()

			def mean(key):
				return float((np.asarray(g[key][:], dtype=float) * w).sum()
					/ w.sum())

			self.energy = mean("energy")
			self.yld = mean("yield")
			self.angle = mean("ion_angle")
			self.hists = {}
			if self.opt == "sheath":
				for key in ["ion_energy_edges", "ion_energy_hist",
						"ion_angle_hist", "sputtering_ion_energy_hist",
						"sputtering_ion_angle_hist"]:
					self.hists[key] = np.asarray(g[key][:], dtype=float)


def moments(s, kmax=5):
	"""
	M_k = int_0^inf x^k exp(-(x - s)^2 / 2) dx for k = 0 .. kmax, by
	quadrature.

	Inputs:
	 - s (float) : The shift, c_s / v_th.

	Returns: list of M_k.
	"""
	x = np.linspace(0.0, s + 12.0, 200001)
	g = np.exp(-0.5 * (x - s) ** 2)
	trapz = getattr(np, "trapezoid", None) or np.trapz
	return [trapz(x ** k * g, x) for k in range(kmax + 1)]


def start_stats(te, ti, num_ions):
	"""
	The flux-weighted mean and standard error of the test ions' starting
	energy.

	Inputs:
	 - te, ti (float)  : Temperatures (eV).
	 - num_ions (int)  : Test ions per point.

	Returns: (mean, standard error) in eV.
	"""
	s = np.sqrt((te + ti) / ti)
	m = moments(s)
	par = 0.5 * ti * m[3] / m[1]
	par_var = (0.5 * ti) ** 2 * (m[5] / m[1] - (m[3] / m[1]) ** 2)
	var = par_var + ti ** 2

	# Weights x: the effective number of samples
	n_eff = num_ions * m[1] ** 2 / (m[0] * m[2])
	return ti + par, np.sqrt(var / n_eff)


def brooks_drop(run):
	"""
	e (phi(start) - phi_wall) / e for Brooks' profile (Sheath::make_params
	and Sheath::potential), with the ions starting 10 rho_s out.

	Inputs:
	 - run (Run) : The run.

	Returns: the drop (V).
	"""
	lam_w = 0.5 * np.log(2.0 * np.pi * run.me / run.mi
		* (1.0 + run.ti / run.te))
	c = abs(np.cos(np.radians(run.theta)))
	ln_cos = np.log(c) if c > 0.0 else -np.inf
	f_d = 1.0 - max(ln_cos, lam_w) / lam_w
	phi_w = lam_w * run.te

	# The Debye part at 10 rho_s is exp(-10 rho_s / (2 lambda_D)), nothing
	phi_start = phi_w * (1.0 - f_d) * np.exp(-START_RHO_S)
	return phi_start - phi_w


def normal_angles(run, n=4000000, seed=1):
	"""
	The flux-weighted angles from the normal the ions hit at with B along
	the normal: the speed across B is kept, the speed along it gains the
	drop.

	Inputs:
	 - run (Run) : The run (theta = 0).

	Returns: (angles in degrees, weights).
	"""
	rng = np.random.default_rng(seed)
	s = np.sqrt((run.te + run.ti) / run.ti)
	x = rng.normal(s, 1.0, n)
	x = x[x > 0.0]
	perp = np.hypot(rng.normal(size=x.size), rng.normal(size=x.size))
	par = np.sqrt(x ** 2 + 2.0 * brooks_drop(run) / run.ti)
	return np.degrees(np.arctan2(perp, par)), x


def hist_stats(counts, centres, n_eff):
	"""
	Mean and standard error from a histogram.

	Returns: (mean, standard error).
	"""
	mean = (counts * centres).sum() / counts.sum()
	var = (counts * (centres - mean) ** 2).sum() / counts.sum()
	return mean, np.sqrt(var / n_eff)


def read_sampler(path):
	"""
	Reads a reg_sheath_angles output: its header values and the counts in
	each 1 degree bin (unweighted).

	Returns: (dict of header values, counts) or None if it isn't there.
	"""
	if not os.path.exists(path):
		return None
	info = {}
	counts = np.zeros(90)
	with open(path) as f:
		for line in f:
			if line.startswith("#"):
				continue
			key, value = line.split()[:2]
			try:
				counts[int(key)] = float(value)
			except ValueError:
				info[key] = float(value)
	return info, counts


def main():
	parser = argparse.ArgumentParser(description=__doc__,
		formatter_class=argparse.RawDescriptionHelpFormatter)
	parser.add_argument("directory", nargs="?", default=".")
	args = parser.parse_args()

	paths = sorted(glob.glob(os.path.join(args.directory,
		"reg_sheath_impacts_*.nc")))
	if not paths:
		print("No reg_sheath_impacts_*.nc in {}".format(args.directory))
		return 1
	runs = {r.name: r for r in map(Run, paths)}
	failed = []

	def check(ok, text):
		print("  {:<4} {}".format("ok" if ok else "FAIL", text))
		if not ok:
			failed.append(text)

	centres = np.arange(90) + 0.5

	# Every sheath run: the table and the histograms
	print("Impact table and histograms")
	for r in runs.values():
		if r.opt != "sheath":
			continue
		a = r.attrs
		check(int(a["ion_impact_lost"]) == 0
			and float(a["ion_impact_max_energy_error"]) < ENERGY_TOL,
			"{:<22} {} lost, energy error {:.1e} of the drop".format(r.name,
			int(a["ion_impact_lost"]), float(a["ion_impact_max_energy_error"])))
		h = r.hists
		ok = all(abs(h[k].sum() / r.ions - 1.0) < HIST_TOL
			for k in ["ion_energy_hist", "ion_angle_hist"])
		ok &= all(abs(h[k].sum() / r.source - 1.0) < HIST_TOL
			for k in ["sputtering_ion_energy_hist",
			"sputtering_ion_angle_hist"])
		check(ok, "{:<22} histograms add up to the ion flux and the "
			"source".format(r.name))

	# Constant yield: the same source either way
	print("\nConstant yield")
	for name, r in runs.items():
		if not name.startswith("const") or r.opt != "sheath":
			continue
		f = runs.get(name + "_formula")
		check(abs(r.yld / CONST_YIELD - 1.0) < SAME_TOL,
			"{:<22} yield {:.10g}".format(name, r.yld))
		if f is not None:
			check(abs(r.source / f.source - 1.0) < SAME_TOL,
				"{:<22} source {:.6e}, formula {:.6e}".format(name, r.source,
				f.source))

	# Linear yield: the mean impact energy
	print("\nYield linear in energy: mean impact energy (eV)")
	print("  {:<4} {:<22} {:>6} {:>6} {:>6} {:>8} {:>8} {:>8} {:>6}".format(
		"", "run", "theta", "Te", "Ti", "drop", "expect", "Flan", "SE"))
	lin_points = []
	for name, r in sorted(runs.items()):
		if not name.startswith("lin"):
			continue
		check(abs(r.yld / (LIN_SLOPE * r.energy) - 1.0) < INTERP_TOL,
			"{:<22} yield {:.6g} = {:g} x energy {:.6g}".format(name, r.yld,
			LIN_SLOPE, r.energy))
		drop = brooks_drop(r)
		e_start, se = start_stats(r.te, r.ti,
			int(r.attrs["ion_impact_num_ions"]))
		expect = drop + e_start
		ok = abs(r.energy - expect) < NUM_SE * se + ENERGY_TOL * drop
		check(ok, "{:<22} {:6.1f} {:6.1f} {:6.2f} {:8.2f} {:8.2f} {:8.2f} "
			"{:6.2f}".format(name, r.theta, r.te, r.ti, drop, expect,
			r.energy, se))
		lin_points.append((name, expect, r.energy, se))

	# B along the normal: the angles
	normal = None
	if "lin_a0" in runs:
		r = runs["lin_a0"]
		ang, w = normal_angles(r)
		mean = np.average(ang, weights=w)
		std = np.sqrt(np.average((ang - mean) ** 2, weights=w))
		n_eff = int(r.attrs["ion_impact_num_ions"]) * (w.sum() ** 2
			/ (w.size * (w ** 2).sum()))
		se = std / np.sqrt(n_eff)
		print("\nB along the normal: flux-weighted mean angle (degrees)")
		check(abs(r.angle - mean) < NUM_SE * se + 0.1,
			"{:<22} expect {:.2f}, Flan {:.2f} (SE {:.2f})".format("lin_a0",
			mean, r.angle, se))
		normal = (r, np.histogram(ang, bins=90, range=(0, 90), weights=w)[0])

	# Borodkina's cases against the sampler at the same plasma
	print("\nBorodkina's cases: flux-weighted mean angle (degrees)")
	bor = {}
	for alpha in [45, 60]:
		r = runs.get("bor_a{}".format(alpha))
		smp = read_sampler(os.path.join(args.directory,
			"reg_sheath_angles_a{}_sampler.txt".format(alpha)))
		if r is None:
			continue
		if smp is None:
			print("  (no reg_sheath_angles_a{}_sampler.txt, not "
				"checked)".format(alpha))
			continue
		info, counts = smp
		n = int(r.attrs["ion_impact_num_ions"])
		_, se_f = hist_stats(r.hists["ion_angle_hist"], centres, 0.8 * n)
		_, se_s = hist_stats(counts, centres, 0.8 * info["hit"])
		se = np.hypot(se_f, se_s)
		check(abs(r.angle - info["mean_angle_flux"]) < NUM_SE * se
			+ GRID_TOL, "bor_a{:<18} sampler {:.2f}, Flan {:.2f} (SE {:.2f}, "
			"grid allowance {:.1f})".format(alpha, info["mean_angle_flux"],
			r.angle, se, GRID_TOL))
		bor[alpha] = (r, info)

	# RustBCA's D on W, printed
	print("\nD on W (RustBCA), formula and sheath: mean energy (eV), angle "
		"(degrees), yield, source (atoms/s)")
	print("  {:>5} {:>5} {:>7} {:>5} {:>9} {:>9} {:>7} {:>5} {:>9} {:>9} "
		"{:>6}".format("Te", "theta", "E", "ang", "Y", "S", "E", "ang", "Y",
		"S", "ratio"))
	dw = []
	for name, r in sorted(runs.items(), key=lambda kv: (kv[1].theta,
			kv[1].te)):
		if not name.startswith("dw") or r.opt != "sheath":
			continue
		f = runs.get(name + "_formula")
		if f is None:
			continue
		ratio = r.source / f.source if f.source > 0 else np.inf
		print("  {:5.0f} {:5.0f} {:7.1f} {:5.1f} {:9.3e} {:9.3e} {:7.1f} "
			"{:5.1f} {:9.3e} {:9.3e} {:6.2f}".format(r.te, r.theta, f.energy,
			f.angle, f.yld, f.source, r.energy, r.angle, r.yld, r.source,
			ratio))
		dw.append((r.theta, r.te, ratio))

	# Plots
	fig, ax = plt.subplots(2, 2, figsize=(11, 8.5))
	a = ax[0, 0]
	for i, (name, expect, flan, se) in enumerate(lin_points):
		a.errorbar(i, flan - expect, yerr=NUM_SE * se, fmt="o", color="C0")
	a.axhline(0.0, color="k", lw=0.8)
	a.set_xticks(range(len(lin_points)))
	a.set_xticklabels([p[0] for p in lin_points], rotation=20)
	a.set_ylabel("Flan - (drop + E_start) (eV)")
	a.set_title("Mean impact energy, linear yield (bars: 4 SE)")

	a = ax[0, 1]
	if normal is not None:
		r, expect = normal
		h = r.hists["ion_angle_hist"]
		a.step(centres, h / h.sum(), where="mid", label="Flan")
		a.plot(centres, expect / expect.sum(), "k--", label="analytic")
		a.legend()
	a.set_xlabel("angle from the normal (degrees)")
	a.set_ylabel("share of the ion flux per degree")
	a.set_title("B along the normal (lin_a0)")

	a = ax[1, 0]
	for i, (alpha, (r, info)) in enumerate(sorted(bor.items())):
		h = r.hists["ion_angle_hist"]
		a.step(centres, h / h.sum(), where="mid", color="C{}".format(i),
			label="Flan, alpha = {}".format(alpha))
		a.axvline(info["mean_angle_flux"], color="C{}".format(i), ls="--",
			label="sampler mean, alpha = {}".format(alpha))
	a.legend(fontsize=8)
	a.set_xlabel("angle from the normal (degrees)")
	a.set_ylabel("share of the ion flux per degree")
	a.set_title("Borodkina's cases")

	a = ax[1, 1]
	for i, theta in enumerate(sorted({d[0] for d in dw})):
		pts = sorted((te, ratio) for th, te, ratio in dw if th == theta)
		a.plot(*zip(*pts), "o-", color="C{}".format(i),
			label="theta = {:g}".format(theta))
	a.axhline(1.0, color="k", lw=0.8)
	a.set_xscale("log")
	a.set_yscale("log")
	tes = sorted({d[1] for d in dw})
	a.set_xticks(tes)
	a.set_xticklabels(["{:g}".format(t) for t in tes])
	a.minorticks_off()
	a.set_xlabel("Te = Ti (eV)")
	a.set_ylabel("source, sheath / formula")
	a.set_title("D on W (RustBCA)")
	a.legend()

	fig.tight_layout()
	out = os.path.join(args.directory, "check_sheath_impacts.png")
	fig.savefig(out, dpi=120)
	print("\nWrote {}".format(out))

	print("\n{} check(s) failed".format(len(failed)) if failed
		else "\nAll checks passed")
	return 1 if failed else 0


if __name__ == "__main__":
	sys.exit(main())
