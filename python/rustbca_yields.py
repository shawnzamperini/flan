"""
Makes a sputtering yield table with RustBCA for Flan.

The yield Y (sputtered target atoms per incident ion) is computed with
RustBCA's sputtering_yield on a grid of incident energies and angles, and
saved to an HDF5 file that Flan reads. The target is a single element,
semi-infinite and homogeneous. The angle is measured from the surface
normal, in degrees.

RustBCA runs every ion with the same seed, and ion i always uses random
stream i, so calling sputtering_yield again with more ions repeats the
first ones. Each point is therefore run twice: a pilot run with n_min ions
to estimate Y, then a run with enough ions to count about target_counts
sputtered atoms (between n_min and n_max ions). If the pilot run sputters
nothing, n_max ions are used. Points that sputter nothing with n_max ions
are stored as Y = 0. Energies are run from highest to lowest, and once an
angle sputters nothing, the lower energies at that angle are below
threshold too, so they're set to 0 without running them (num_samples = 0).
The error is the counting error sqrt(counts) / N.

The file holds:

 - /energy [nE] : Incident energy (eV).
 - /angle [nA] : Incident angle from the surface normal (degrees).
 - /yield [nE, nA] : Sputtering yield (atoms/ion).
 - /yield_err [nE, nA] : Counting error of the yield (atoms/ion).
 - /num_samples [nE, nA] : Number of ions run for each point.
 - /projectile/{Z, mass_amu, Es, Ec} : Projectile parameters.
 - /target/{Z, mass_amu, n, Es, Eb, Ec} : Target parameters. n is in m^-3,
   the energies in eV.

Provenance (RustBCA version and commit, seed, date, command line and the
physics models RustBCA uses) is stored as attributes of the root group.

Usage from the command line (the flan conda environment has RustBCA):

	python rustbca_yields.py D W
	python rustbca_yields.py D W --quick --plot
	python rustbca_yields.py --plot-only ../data/sputtering/D_on_W.h5

The default output is ../data/sputtering/{projectile}_on_{target}.h5,
relative to this script.
"""

import argparse
import datetime
import importlib.metadata
import os
import pathlib
import subprocess
import sys
import time

import h5py
import numpy as np


# Projectile and target parameters. Z is the atomic number, m the mass
# (amu), Es the surface binding energy (eV), Ec the cutoff energy below
# which RustBCA stops following a particle (eV), Eb the bulk binding energy
# (eV) and n the number density (m^-3). Only targets need Eb and n. Es for
# W is its sublimation energy, the same binding energy Flan uses for the
# Thompson distribution of sputtered atoms. The rest follow RustBCA's
# scripts/materials.py.
MATERIALS = {
	"D": {"Z": 1.0, "m": 2.014, "Es": 1.5, "Ec": 0.1},
	"W": {"Z": 74.0, "m": 183.84, "Es": 8.79, "Ec": 3.0, "Eb": 0.0,
		"n": 6.306e28},
}

# Physics models sputtering_yield uses, from RustBCA's
# Options::default_options at the pinned commit. They can't be changed
# through sputtering_yield, so they're only recorded.
RUSTBCA_MODELS = {
	"interaction_potential": "KR_C",
	"electronic_stopping_mode": "LOW_ENERGY_NONLOCAL",
	"mean_free_path_model": "LIQUID",
	"scattering_integral": "MENDENHALL_WELLER",
	"surface_binding_model": "AVERAGE",
	"bulk_binding_model": "INDIVIDUAL",
}

# Angles used by default (degrees from the normal). Finer near grazing,
# where the yield changes fastest.
DEFAULT_ANGLES = list(range(0, 90, 5)) + [87, 89]


def rustbca_commit():
	"""
	Returns the commit of the RustBCA clone in $FLANSOFT/RustBCA, or
	"unknown" if it can't be found.
	"""
	flansoft = os.environ.get("FLANSOFT",
		os.path.join(os.path.expanduser("~"), "flansoft"))
	try:
		return subprocess.run(["git", "-C",
			os.path.join(flansoft, "RustBCA"), "rev-parse", "HEAD"],
			capture_output=True, text=True, check=True).stdout.strip()
	except (OSError, subprocess.CalledProcessError):
		return "unknown"


def yield_point(sputtering_yield, ion, target, energy, angle,
	target_counts, n_min, n_max):
	"""
	Computes the yield at one energy and angle.

	Inputs
	 - sputtering_yield (function) : libRustBCA.sputtering_yield.
	 - ion (dict)           : Projectile parameters (Z, m, Es, Ec).
	 - target (dict)        : Target parameters (Z, m, Es, Ec, Eb, n).
	 - energy (float)       : Incident energy (eV).
	 - angle (float)        : Incident angle from the normal (degrees).
	 - target_counts (int)  : Number of sputtered atoms to aim for.
	 - n_min (int)          : Ions in the pilot run, and the fewest used.
	 - n_max (int)          : Most ions used.

	Returns the yield, its counting error and the number of ions run.
	"""
	y = sputtering_yield(ion, target, energy, angle, n_min)
	if y > 0.0:
		n = int(np.clip(np.ceil(target_counts / y), n_min, n_max))
	else:
		n = n_max
	if n != n_min:
		y = sputtering_yield(ion, target, energy, angle, n)
	err = np.sqrt(y * n) / n
	return y, err, n


def make_table(projectile, target, energies, angles, target_counts,
	n_min, n_max, seed):
	"""
	Computes the yield table.

	Inputs
	 - projectile (dict)   : Projectile parameters (Z, m, Es, Ec).
	 - target (dict)       : Target parameters (Z, m, Es, Ec, Eb, n).
	 - energies (array)    : Incident energies (eV).
	 - angles (array)      : Incident angles from the normal (degrees).
	 - target_counts (int) : Number of sputtered atoms to aim for per point.
	 - n_min (int)         : Ions in each pilot run.
	 - n_max (int)         : Most ions per point.
	 - seed (int)          : RustBCA seed (LIBRUSTBCA_SEED).

	Returns the yield, its error and the number of ions run, each with
	shape (len(energies), len(angles)).
	"""
	# RustBCA reads the seed from the environment on every call
	os.environ["LIBRUSTBCA_SEED"] = str(seed)
	from libRustBCA import sputtering_yield

	shape = (len(energies), len(angles))
	yields = np.zeros(shape)
	errs = np.zeros(shape)
	nums = np.zeros(shape, dtype=np.int64)

	# Energies are run from highest to lowest. Once an angle sputters
	# nothing with n_max ions, it's below threshold, and so are the lower
	# energies, which are left at 0 without running them.
	below = np.zeros(len(angles), dtype=bool)
	start = time.time()
	for i in np.argsort(energies)[::-1]:
		for j, angle in enumerate(angles):
			if below[j]:
				continue
			yields[i, j], errs[i, j], nums[i, j] = yield_point(
				sputtering_yield, projectile, target, float(energies[i]),
				float(angle), target_counts, n_min, n_max)
			below[j] = yields[i, j] == 0.0

		print("E = {:9.2f} eV: max Y = {:.3e}, {:6.0f} s elapsed".format(
			energies[i], yields[i].max(), time.time() - start), flush=True)
		if below.all():
			print("Every angle is below threshold, skipping lower energies")
			break

	return yields, errs, nums


def save_table(path, energies, angles, yields, errs, nums, projectile,
	target, projectile_name, target_name, seed):
	"""
	Saves the yield table and its provenance to an HDF5 file.

	Inputs
	 - path (str)            : Output file.
	 - energies (array)      : Incident energies (eV).
	 - angles (array)        : Incident angles from the normal (degrees).
	 - yields (array)        : Yields (atoms/ion), shape (nE, nA).
	 - errs (array)          : Errors of the yields, shape (nE, nA).
	 - nums (array)          : Ions run per point, shape (nE, nA).
	 - projectile (dict)     : Projectile parameters.
	 - target (dict)         : Target parameters.
	 - projectile_name (str) : Projectile name, e.g. "D".
	 - target_name (str)     : Target name, e.g. "W".
	 - seed (int)            : RustBCA seed used.
	"""
	pathlib.Path(path).parent.mkdir(parents=True, exist_ok=True)
	with h5py.File(path, "w") as f:
		f.create_dataset("energy", data=energies)
		f.create_dataset("angle", data=angles)
		f.create_dataset("yield", data=yields)
		f.create_dataset("yield_err", data=errs)
		f.create_dataset("num_samples", data=nums)

		proj = f.create_group("projectile")
		proj.create_dataset("Z", data=projectile["Z"])
		proj.create_dataset("mass_amu", data=projectile["m"])
		proj.create_dataset("Es", data=projectile["Es"])
		proj.create_dataset("Ec", data=projectile["Ec"])

		targ = f.create_group("target")
		targ.create_dataset("Z", data=target["Z"])
		targ.create_dataset("mass_amu", data=target["m"])
		targ.create_dataset("n", data=target["n"])
		targ.create_dataset("Es", data=target["Es"])
		targ.create_dataset("Eb", data=target["Eb"])
		targ.create_dataset("Ec", data=target["Ec"])

		f["energy"].attrs["units"] = "eV"
		f["angle"].attrs["units"] = "degrees from the surface normal"
		f["yield"].attrs["units"] = "sputtered atoms per incident ion"
		f["target/n"].attrs["units"] = "m^-3"

		f.attrs["projectile"] = projectile_name
		f.attrs["target"] = target_name
		f.attrs["rustbca_version"] = importlib.metadata.version("RustBCA")
		f.attrs["rustbca_commit"] = rustbca_commit()
		f.attrs["seed"] = seed
		f.attrs["created"] = datetime.datetime.now().isoformat(
			timespec="seconds")
		f.attrs["command"] = " ".join(sys.argv)
		for key, val in RUSTBCA_MODELS.items():
			f.attrs[key] = val


def plot_table(path, save_prefix=None):
	"""
	Plots a yield table: Y(E) at a few angles, and Y(angle) at a few
	energies.

	Inputs
	 - path (str)        : Yield table file.
	 - save_prefix (str) : If given, the plots are saved as
	                       {save_prefix}_energy.png and
	                       {save_prefix}_angle.png instead of shown.
	"""
	import matplotlib.pyplot as plt

	with h5py.File(path, "r") as f:
		energies = f["energy"][:]
		angles = f["angle"][:]
		yields = f["yield"][:]
		errs = f["yield_err"][:]
		title = "{} on {}".format(f.attrs["projectile"], f.attrs["target"])

	# Table points nearest to the values asked for, without repeats (a
	# coarse table can have the same nearest point for two of them)
	def nearest(grid, wanted):
		return list(dict.fromkeys(np.argmin(np.abs(grid - w))
			for w in wanted))

	# Y(E) at a few angles. Zero yields can't go on a log axis.
	fig1, ax1 = plt.subplots(figsize=(6, 4.5))
	for j in nearest(angles, [0, 30, 60, 80, 89]):
		mask = yields[:, j] > 0
		ax1.errorbar(energies[mask], yields[mask, j], yerr=errs[mask, j],
			marker=".", capsize=2, label="{:g}°".format(angles[j]))
	ax1.set_xscale("log")
	ax1.set_yscale("log")
	ax1.set_xlabel("Incident energy (eV)")
	ax1.set_ylabel("Sputtering yield (atoms/ion)")
	ax1.set_title(title)
	ax1.legend(title="Angle from normal")
	ax1.grid(which="both", alpha=0.3)
	fig1.tight_layout()

	# Y(angle) at a few energies
	fig2, ax2 = plt.subplots(figsize=(6, 4.5))
	for i in nearest(energies, [300, 500, 1000, 3000, 10000]):
		ax2.errorbar(angles, yields[i], yerr=errs[i], marker=".",
			capsize=2, label="{:.0f} eV".format(energies[i]))
	ax2.set_yscale("log")
	ax2.set_xlabel("Angle from normal (degrees)")
	ax2.set_ylabel("Sputtering yield (atoms/ion)")
	ax2.set_title(title)
	ax2.legend(title="Energy")
	ax2.grid(which="both", alpha=0.3)
	fig2.tight_layout()

	if save_prefix is not None:
		fig1.savefig(save_prefix + "_energy.png", dpi=150)
		fig2.savefig(save_prefix + "_angle.png", dpi=150)
		print("Saved {}_energy.png and {}_angle.png".format(save_prefix,
			save_prefix))
	else:
		plt.show()


def main():
	parser = argparse.ArgumentParser(description="Make a sputtering yield "
		"table with RustBCA for Flan.")
	parser.add_argument("projectile", nargs="?", default="D",
		choices=[k for k, v in MATERIALS.items()],
		help="Incident ion (default D)")
	parser.add_argument("target", nargs="?", default="W",
		choices=[k for k, v in MATERIALS.items() if "n" in v],
		help="Target material (default W)")
	parser.add_argument("--out", help="Output file (default "
		"../data/sputtering/{projectile}_on_{target}.h5 relative to this "
		"script)")
	parser.add_argument("--emin", type=float, default=10.0,
		help="Lowest energy, eV (default 10)")
	parser.add_argument("--emax", type=float, default=1e4,
		help="Highest energy, eV (default 1e4)")
	parser.add_argument("--ne", type=int, default=41,
		help="Number of log-spaced energies (default 41)")
	parser.add_argument("--angles", type=float, nargs="+",
		default=DEFAULT_ANGLES, help="Angles from the normal, degrees "
		"(default 0 to 85 every 5, 87 and 89)")
	parser.add_argument("--target-counts", type=int, default=2000,
		help="Sputtered atoms to aim for per point (default 2000, about "
		"2%% error)")
	parser.add_argument("--n-min", type=int, default=10000,
		help="Ions in each pilot run, and the fewest per point (default "
		"1e4)")
	parser.add_argument("--n-max", type=int, default=2000000,
		help="Most ions per point (default 2e6)")
	parser.add_argument("--seed", type=int, default=0,
		help="RustBCA seed (default 0)")
	parser.add_argument("--target-Es", type=float,
		help="Override the target's surface binding energy, eV")
	parser.add_argument("--quick", action="store_true",
		help="Coarse grid and few ions, for testing")
	parser.add_argument("--plot", action="store_true",
		help="Plot the table after making it")
	parser.add_argument("--plot-only", metavar="FILE",
		help="Only plot an existing table")
	parser.add_argument("--save-plots", metavar="PREFIX",
		help="Save the plots as PREFIX_energy.png and PREFIX_angle.png "
		"instead of showing them")
	args = parser.parse_args()

	if args.plot_only:
		plot_table(args.plot_only, args.save_plots)
		return

	projectile = dict(MATERIALS[args.projectile])
	target = dict(MATERIALS[args.target])
	if args.target_Es is not None:
		target["Es"] = args.target_Es

	energies = np.geomspace(args.emin, args.emax, args.ne)
	angles = np.array(args.angles, dtype=float)
	target_counts, n_min, n_max = args.target_counts, args.n_min, args.n_max
	if args.quick:
		energies = np.geomspace(args.emin, args.emax, 9)
		angles = np.array([0.0, 30.0, 60.0, 85.0])
		target_counts, n_min, n_max = 200, 2000, 20000
	if np.any(angles < 0) or np.any(angles > 90):
		parser.error("Angles must be between 0 and 90 degrees.")

	out = args.out
	if out is None:
		out = str(pathlib.Path(__file__).resolve().parent.parent / "data"
			/ "sputtering" / "{}_on_{}.h5".format(args.projectile,
			args.target))

	print("{} on {}: {} energies x {} angles, saving to {}".format(
		args.projectile, args.target, len(energies), len(angles), out))
	yields, errs, nums = make_table(projectile, target, energies, angles,
		target_counts, n_min, n_max, args.seed)
	save_table(out, energies, angles, yields, errs, nums, projectile,
		target, args.projectile, args.target, args.seed)
	print("Saved {}".format(out))

	if args.plot:
		plot_table(out, args.save_plots)


if __name__ == "__main__":
	main()
