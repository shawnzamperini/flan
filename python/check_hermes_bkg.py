"""
Checks a HERMES background the way Flan used it, from a Flan output file.

python/check_hermes_conventions.py checks the HERMES files themselves. This
script reads what Flan made of them (the background saved in the output
file), so it also checks the reader's signs, units and geometry:

1. Flow into the targets. At each target (a y face that is a wall), the
   parallel ion flow u_par = u . b, taken toward the target, is compared
   with the sound speed c_s = sqrt(e (Te + Ti) / m_i). The sheath condition
   of the run makes it about sonic or faster, so a median well under 1, or
   flow away from the target, means a sign or unit is wrong.
2. Grad-B drift direction. The ion grad-B drift, along B x grad B, is
   found at the outer midplane just outside the separatrix and compared
   with the X-point's position. With the conventional "forward" field the
   ion grad-B drift points toward the X-point; reversing the field
   reverses it. Give --expect to make it a pass/fail check.
3. Impurity content (with --content). The impurity content of each zone,
   from the density times the cell volumes, as a fraction of the total.
   Comparing two runs (e.g., the two field directions, or
   hermes_field_aligned on and off) shows where the impurities go.

The zones are numbered as Flan's HERMES reader makes them: 1 inner PFR,
2 core, 3 outer PFR, 4 inner leg, 5 main SOL, 6 outer leg.

Usage:

	python check_hermes_bkg.py hermes_example.nc
	python check_hermes_bkg.py hermes_example.nc --expect toward --content

The exit code is 1 if a check fails, 0 otherwise.
"""

import argparse
import sys

import netCDF4
import numpy as np

QE = 1.602176634e-19
MP = 1.67262192369e-27

# Flan's zone face codes (see the zones group)
FACE_WALL = -4
ZONE_NAMES = {1: "inner PFR", 2: "core", 3: "outer PFR", 4: "inner leg",
	5: "main SOL", 6: "outer leg"}


def vec(group, name):
	"""Cartesian vector (name_X, name_Y, name_Z) or (nameX, ...) as an
	array with the components in the last axis."""
	keys = [k for k in (name + "_X", name + "X") if k in group.variables]
	if not keys:
		raise KeyError(name)
	base = keys[0][:-1]
	return np.stack([np.asarray(group[base + c][:]) for c in "XYZ"], axis=-1)


def check_targets(nc, mi):
	"""
	Prints the Mach number of the parallel ion flow into each target and
	returns True if every target's median is at least 0.8 toward it.

	Inputs
	 - nc (netCDF4.Dataset) : Open Flan output file.
	 - mi (float)           : Main ion mass (kg).
	"""
	faces = np.asarray(nc["zones"]["faces"][:])
	ok = True
	print("Parallel ion flow into the targets (u_par toward the target / c_s):")
	for zone in range(1, faces.shape[1] + 1):
		for face, end in ((2, 0), (3, -1)):
			if faces[face, zone - 1] != FACE_WALL:
				continue
			g = nc["zone{}".format(zone)]
			bkg, geo = g["background"], g["geometry"]
			B = vec(bkg, "B")[:, :, end]
			U = vec(bkg, "Ui")[:, :, end]
			grad_y = np.stack([np.asarray(geo["dyd" + c][:, end])
				for c in "XYZ"], axis=-1)
			Bm = np.linalg.norm(B, axis=-1)
			u_par = np.sum(U * B, axis=-1) / Bm

			# Which way along B is +y, and which way the target is
			along_y = np.sign(np.sum(B * grad_y[None], axis=-1))
			toward = 1.0 if face == 3 else -1.0
			te = np.asarray(bkg["Te"][:, :, end])
			ti = np.asarray(bkg["Ti"][:, :, end])
			cs = np.sqrt(QE * np.maximum(te + ti, 0) / mi)
			good = cs > 0
			mach = (u_par * along_y * toward)[good] / cs[good]
			med = np.median(mach)
			print("  zone {} ({}), {} target: median {:+.2f}, 10th "
				"percentile {:+.2f}, {:.0f}% of cells >= 0.9".format(zone,
				ZONE_NAMES.get(zone, "?"), "ymax" if face == 3 else "ymin",
				med, np.percentile(mach, 10), 100 * np.mean(mach >= 0.9)))
			ok = ok and med >= 0.8
	return ok


def check_gradb(nc, expect):
	"""
	Prints which way the ion grad-B drift points at the outer midplane,
	relative to the X-point, and returns False if it isn't expect
	("toward" or "away"; None to only print it).
	"""
	sol = nc["zone5"]
	bkg, geo = sol["background"], sol["geometry"]
	R = np.hypot(np.asarray(geo["X"][0]), np.asarray(geo["Y"][0]))
	Zc = np.asarray(geo["Z"][0])
	j = int(np.argmax(R.mean(axis=-1)))
	B = vec(bkg, "B")[:, 0, j].mean(axis=0)
	gB = vec(bkg, "dBd")[:, 0, j].mean(axis=0)
	drift = np.cross(B, gB)
	drift_Z = drift[:, 2].mean()

	core = nc["zone2"]["geometry"]
	z_x = float(core["Z_corners"][-1, 0])
	z_mid = float(Zc[j].mean())
	toward = np.sign(drift_Z) == np.sign(z_x - z_mid)
	print("Ion grad-B drift at the outer midplane (R = {:.3f} m): points {} "
		"({} the X-point at Z = {:.3f} m)".format(R[j].mean(),
		"down" if drift_Z < 0 else "up", "toward" if toward else "away from",
		z_x))
	if expect is None:
		return True
	ok = (expect == "toward") == toward
	print("  expected {}: {}".format(expect, "pass" if ok else "FAIL"))
	return ok


def print_content(nc):
	"""Prints each zone's share of the impurities, averaged over time."""
	nzones = len(nc["zones"]["nx"][:])
	content = []
	for zone in range(1, nzones + 1):
		g = nc["zone{}".format(zone)]
		nz = np.asarray(g["output"]["nz"][:]).mean(axis=0)
		J = np.asarray(g["geometry"]["J"][:])
		zc = np.asarray(g["geometry"]["z"][:])
		dz = zc[1] - zc[0] if zc.size > 1 else 1.0
		content.append(np.sum(nz * J) * dz)
	total = sum(content)
	print("Impurity content by zone (fraction of the total):")
	for zone, c in enumerate(content, start=1):
		print("  zone {} ({}): {:.3f}".format(zone, ZONE_NAMES.get(zone, "?"),
			c / total if total > 0 else 0.0))


def main():
	parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0],
		formatter_class=argparse.RawDescriptionHelpFormatter)
	parser.add_argument("path", help="Flan output file of a HERMES run")
	parser.add_argument("--expect", choices=["toward", "away"],
		help="Expected direction of the ion grad-B drift relative to the "
		"X-point: toward for the forward field, away for the reversed")
	parser.add_argument("--content", action="store_true",
		help="Also print the impurity content of each zone")
	args = parser.parse_args()

	with netCDF4.Dataset(args.path) as nc:
		if nc["input"]["bkg_source"][0] != "hermes":
			print("Warning: bkg_source isn't \"hermes\".")
		mi = float(nc["input"]["hermes_ion_mass_amu"][0]) * MP
		ok = check_targets(nc, mi)
		ok = check_gradb(nc, args.expect) and ok
		if args.content:
			print_content(nc)
	print("All checks passed" if ok else "Some checks FAILED")
	return 0 if ok else 1


if __name__ == "__main__":
	sys.exit(main())
