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

The yield rises steeply near threshold, too steeply to interpolate
accurately on a coarse grid, so after the main grid is run, energies about
5% apart are added from the highest energy that sputtered nothing up to 3
times the lowest energy that sputtered something (--fine-step,
--fine-factor, or --no-refine to turn this off).

Every point is independent and RustBCA is seeded the same way each time,
so a point comes out the same in every table made with the same settings.
--reuse copies the points of an existing table made with the same
settings instead of running them again, e.g., to add energies to a table.

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


def point_key(energy, angle):
	"""
	Returns a key for a point of a table, the same for an energy and angle
	that are equal to 9 significant figures, so points from a saved table
	can be matched to the ones being computed.
	"""
	return ("{:.9g}".format(energy), "{:.9g}".format(angle))


def run_energies(sputtering_yield, projectile, target, energies, angles,
	target_counts, n_min, n_max, reuse, start):
	"""
	Computes the yield at each of some energies, for every angle.

	Energies are run from highest to lowest. Once an angle sputters nothing
	with n_max ions, it's below threshold, and so are the lower energies,
	which are left at 0 without running them (num_samples = 0).

	Inputs
	 - sputtering_yield (function) : libRustBCA.sputtering_yield.
	 - projectile (dict)   : Projectile parameters (Z, m, Es, Ec).
	 - target (dict)       : Target parameters (Z, m, Es, Ec, Eb, n).
	 - energies (array)    : Incident energies (eV), in any order.
	 - angles (array)      : Incident angles from the normal (degrees).
	 - target_counts (int) : Number of sputtered atoms to aim for per point.
	 - n_min (int)         : Ions in each pilot run.
	 - n_max (int)         : Most ions per point.
	 - reuse (dict)        : Points already computed with the same
	                         settings, from point_key to (yield, error,
	                         ions), used instead of running them again.
	 - start (float)       : Time the table was started, for progress.

	Returns a dictionary from each energy to its yields, errors and ions
	run, each an array over the angles.
	"""
	rows = {}
	below = np.zeros(len(angles), dtype=bool)
	for energy in sorted(energies, reverse=True):
		y = np.zeros(len(angles))
		err = np.zeros(len(angles))
		num = np.zeros(len(angles), dtype=np.int64)
		reused = 0
		for j, angle in enumerate(angles):
			if below[j]:
				continue
			key = point_key(energy, angle)
			if key in reuse:
				y[j], err[j], num[j] = reuse[key]
				reused += 1
			else:
				y[j], err[j], num[j] = yield_point(sputtering_yield,
					projectile, target, float(energy), float(angle),
					target_counts, n_min, n_max)
			below[j] = y[j] == 0.0
		rows[energy] = (y, err, num)

		print("E = {:9.2f} eV: max Y = {:.3e}, {:6.0f} s elapsed{}".format(
			energy, y.max(), time.time() - start,
			", {} points reused".format(reused) if reused else ""),
			flush=True)
		if below.all():
			print("Every angle is below threshold, skipping lower energies")
			break
	return rows


def threshold_energies(rows, emax, fine_step, fine_factor):
	"""
	Picks extra energies near threshold, where the yield changes fastest.

	The threshold region runs from the highest energy where every angle
	sputtered nothing (or the lowest energy, if none did) up to fine_factor
	times the lowest energy where any angle sputtered something, but no
	higher than emax. Energies about fine_step apart (as a fraction) are
	put across it, leaving out those within fine_step / 4 of an energy
	that's already in rows.

	Inputs
	 - rows (dict)         : Energies computed so far, from run_energies.
	                         Energies missing from it were below threshold.
	 - emax (float)        : Highest energy of the table (eV).
	 - fine_step (float)   : Fractional spacing of the extra energies.
	 - fine_factor (float) : How far above the lowest energy that
	                         sputtered something the region goes.

	Returns the extra energies (eV), which may be empty.
	"""
	energies = np.array(sorted(rows))
	nonzero = [e for e in energies if rows[e][0].max() > 0.0]
	if not nonzero:
		return np.array([])
	first = min(nonzero)
	zeros = [e for e in energies if e < first]
	e_lo = max(zeros) if zeros else first
	e_hi = min(fine_factor * first, emax)
	if e_hi <= e_lo:
		return np.array([])

	n = int(np.ceil(np.log(e_hi / e_lo) / np.log(1.0 + fine_step)))
	fine = np.geomspace(e_lo, e_hi, n + 1)
	keep = [e for e in fine
		if np.min(np.abs(energies / e - 1.0)) > fine_step / 4.0]
	return np.array(keep)


def load_reuse(path, projectile, target, seed, target_counts, n_min,
	n_max):
	"""
	Loads the points of an existing table so they aren't run again.

	The table must have been made with the same projectile and target
	parameters, seed, target_counts, n_min and n_max, or its points would
	be different from the ones being computed, so the script stops if any
	differ. Points that weren't run (num_samples = 0) aren't loaded.

	Inputs
	 - path (str)          : Existing table.
	 - projectile (dict)   : Projectile parameters being used.
	 - target (dict)       : Target parameters being used.
	 - seed (int)          : RustBCA seed being used.
	 - target_counts (int) : Sputtered atoms aimed for per point.
	 - n_min (int)         : Ions in each pilot run.
	 - n_max (int)         : Most ions per point.

	Returns a dictionary from point_key to (yield, error, ions).
	"""
	with h5py.File(path, "r") as f:
		wanted = {"projectile/Z": projectile["Z"],
			"projectile/mass_amu": projectile["m"],
			"projectile/Es": projectile["Es"],
			"projectile/Ec": projectile["Ec"],
			"target/Z": target["Z"], "target/mass_amu": target["m"],
			"target/n": target["n"], "target/Es": target["Es"],
			"target/Eb": target["Eb"], "target/Ec": target["Ec"]}
		for name, val in wanted.items():
			if f[name][()] != val:
				sys.exit("Can't reuse {}: {} is {}, not {}.".format(path,
					name, f[name][()], val))
		settings = {"seed": seed, "target_counts": target_counts,
			"n_min": n_min, "n_max": n_max}
		for name, val in settings.items():
			if name not in f.attrs or f.attrs[name] != val:
				sys.exit("Can't reuse {}: its {} is {}, not {}.".format(
					path, name, f.attrs.get(name, "not recorded"), val))

		energies = f["energy"][:]
		angles = f["angle"][:]
		yields = f["yield"][:]
		errs = f["yield_err"][:]
		nums = f["num_samples"][:]

	reuse = {}
	for i, energy in enumerate(energies):
		for j, angle in enumerate(angles):
			if nums[i, j] > 0:
				reuse[point_key(energy, angle)] = (yields[i, j], errs[i, j],
					nums[i, j])
	print("Reusing {} points from {}".format(len(reuse), path))
	return reuse


def make_table(projectile, target, energies, angles, target_counts,
	n_min, n_max, seed, refine=True, fine_step=0.05, fine_factor=3.0,
	reuse=None):
	"""
	Computes the yield table.

	The energies given are run first. Then, if refine is on, extra
	energies are added near threshold (see threshold_energies) and run
	too, since the yield rises steeply there and a coarse grid can't be
	interpolated accurately.

	Inputs
	 - projectile (dict)   : Projectile parameters (Z, m, Es, Ec).
	 - target (dict)       : Target parameters (Z, m, Es, Ec, Eb, n).
	 - energies (array)    : Incident energies (eV).
	 - angles (array)      : Incident angles from the normal (degrees).
	 - target_counts (int) : Number of sputtered atoms to aim for per point.
	 - n_min (int)         : Ions in each pilot run.
	 - n_max (int)         : Most ions per point.
	 - seed (int)          : RustBCA seed (LIBRUSTBCA_SEED).
	 - refine (bool)       : Add energies near threshold.
	 - fine_step (float)   : Fractional spacing of the added energies.
	 - fine_factor (float) : How far above the threshold to add energies,
	                         as a multiple of the lowest energy that
	                         sputtered something.
	 - reuse (dict)        : Points already computed, from load_reuse, or
	                         None.

	Returns the energies (sorted), and the yield, its error and the number
	of ions run, each with shape (number of energies, len(angles)).
	"""
	# RustBCA reads the seed from the environment on every call
	os.environ["LIBRUSTBCA_SEED"] = str(seed)
	from libRustBCA import sputtering_yield

	reuse = reuse or {}
	start = time.time()
	rows = run_energies(sputtering_yield, projectile, target, energies,
		angles, target_counts, n_min, n_max, reuse, start)

	all_energies = list(energies)
	if refine:
		fine = threshold_energies(rows, max(energies), fine_step,
			fine_factor)
		if len(fine):
			print("Adding {} energies near threshold, {:.1f} to {:.1f} "
				"eV".format(len(fine), fine.min(), fine.max()), flush=True)
			rows.update(run_energies(sputtering_yield, projectile, target,
				fine, angles, target_counts, n_min, n_max, reuse, start))
			all_energies += list(fine)

	all_energies = np.array(sorted(all_energies))
	shape = (len(all_energies), len(angles))
	yields = np.zeros(shape)
	errs = np.zeros(shape)
	nums = np.zeros(shape, dtype=np.int64)
	for i, energy in enumerate(all_energies):
		if energy in rows:
			yields[i], errs[i], nums[i] = rows[energy]

	return all_energies, yields, errs, nums


def save_table(path, energies, angles, yields, errs, nums, projectile,
	target, projectile_name, target_name, settings):
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
	 - settings (dict)       : How the table was made (seed,
	                           target_counts, n_min, n_max, refine,
	                           fine_step, fine_factor), saved as root
	                           attributes. load_reuse checks the first
	                           four.
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
		for key, val in settings.items():
			f.attrs[key] = val
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
	parser.add_argument("--no-refine", action="store_true",
		help="Don't add energies near threshold")
	parser.add_argument("--fine-step", type=float, default=0.05,
		help="Fractional spacing of the energies added near threshold "
		"(default 0.05)")
	parser.add_argument("--fine-factor", type=float, default=3.0,
		help="Energies are added near threshold up to this times the "
		"lowest energy that sputters anything (default 3)")
	parser.add_argument("--reuse", metavar="FILE",
		help="Copy points from an existing table made with the same "
		"materials, seed, --target-counts, --n-min and --n-max, instead "
		"of running them again")
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

	reuse = None
	if args.reuse:
		reuse = load_reuse(args.reuse, projectile, target, args.seed,
			target_counts, n_min, n_max)

	print("{} on {}: {} energies x {} angles, saving to {}".format(
		args.projectile, args.target, len(energies), len(angles), out))
	energies, yields, errs, nums = make_table(projectile, target, energies,
		angles, target_counts, n_min, n_max, args.seed,
		refine=not args.no_refine, fine_step=args.fine_step,
		fine_factor=args.fine_factor, reuse=reuse)
	settings = {"seed": args.seed, "target_counts": target_counts,
		"n_min": n_min, "n_max": n_max, "refine": int(not args.no_refine),
		"fine_step": args.fine_step, "fine_factor": args.fine_factor}
	save_table(out, energies, angles, yields, errs, nums, projectile,
		target, args.projectile, args.target, settings)
	print("Saved {}".format(out))

	if args.plot:
		plot_table(out, args.save_plots)


if __name__ == "__main__":
	main()
