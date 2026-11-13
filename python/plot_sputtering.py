"""
Summarizes and plots a Flan run with imp_source_opt = "sputter": where the
wall is eroded and where the sputtered atoms come back, how quickly and in
what charge state, and where the impurity ends up.

It reads the run's netCDF file (the sputtering group, and the zones for the
impurity inventory) and prints:

 - the sputtered source, the share redeposited, and what happened to the
   particles (fates);
 - the mean time a particle was followed, and the impurity inventory per
   zone and away from the divertor, checked against source x time;

and plots:

 - the gross erosion and net erosion per unit area along the wall in the
   (R, Z) plane, summed over the toroidal cells;
 - how long after being sputtered the redeposited atoms landed, with the
   gyroperiod of a singly charged ion at the strongest source marked;
 - their charge when they landed.

Usage from the command line:

	python plot_sputtering.py s3x_example.nc
	python plot_sputtering.py s3x_example.nc --mesh /path/to/mesh_S3X.h5 --save run1

or from Python:

	from plot_sputtering import load_run, print_summary, plot_all
	run = load_run("s3x_example.nc")
	print_summary(run)
	plot_all(run)

The inventory away from the divertor assumes a single null: it's the part
above the X-point for a lower single null, whose height is taken as the
lowest corner of the core zones (the zones with a core face), or below it
for an upper single null (--upper), using their highest corner.
"""

import argparse

import matplotlib.pyplot as plt
import netCDF4
import numpy as np
from matplotlib.colors import LogNorm, SymLogNorm

# Elementary charge (C) and atomic mass unit (kg)
EV = 1.602176634e-19
AMU = 1.66053906660e-27


def load_run(nc_path):
	"""
	Loads what the summary and plots need from a Flan output file.

	Inputs
	 - nc_path (str) : Flan output file from a run with
	                   imp_source_opt = "sputter".

	Returns a dictionary with the sputtering group's variables and
	attributes, the run's time step, number of particles and impurity
	mass, and per zone the impurity density averaged over the run, the
	cell volumes, the cell-center Z, the number of time steps recorded and
	|B| in the first frame.

	Each frame's nz is the density near that frame's time, the average
	over the times nearer to it than to any other frame. The run average
	weights each frame by the length of that window.
	"""
	run = {}
	with netCDF4.Dataset(nc_path) as nc:
		if "sputtering" not in nc.groups:
			raise ValueError("{} has no sputtering group. Run with "
				"imp_source_opt = \"sputter\".".format(nc_path))
		g = nc["sputtering"]
		for name in g.variables:
			run[name] = g[name][:]
		run["attrs"] = {a: g.getncattr(a) for a in g.ncattrs()}

		inp = nc["input"]
		run["dt"] = float(np.ravel(inp["imp_time_step"][:])[0])
		run["imp_num"] = int(np.ravel(inp["imp_num"][:])[0])
		run["mass_amu"] = float(np.ravel(inp["imp_mass_amu"][:])[0])

		# Zones with a core face (code -2), for the X-point height
		faces = nc["zones"]["faces"][:]
		nzones = faces.shape[1]
		run["core_zones"] = [z + 1 for z in range(nzones)
			if (faces[:, z] == -2).any()]

		# Each frame's window: half way to its neighbors, and to the ends
		times = np.asarray(nc["geometry"]["time"][:], dtype=float)
		if len(times) > 1:
			mids = 0.5 * (times[1:] + times[:-1])
			lo = np.concatenate([[times[0]], mids])
			hi = np.concatenate([mids, [times[-1]]])
			window = hi - lo
		else:
			window = np.ones(1)
		w = window / window.sum()

		run["zones"] = {}
		for z in range(1, nzones + 1):
			zg = nc["zone{}".format(z)]
			geo = zg["geometry"]
			gx, gy, gz = geo["grid_x"][:], geo["grid_y"][:], geo["grid_z"][:]
			vol = (geo["J"][:] * np.diff(gx)[:, None, None]
				* np.diff(gy)[None, :, None] * np.diff(gz)[None, None, :])

			# Cell-center Z from the corners (first toroidal plane)
			Zc = geo["Z_corners"][:]
			Z_cen = 0.25 * (Zc[:-1, :-1] + Zc[1:, :-1] + Zc[:-1, 1:]
				+ Zc[1:, 1:])

			bkg = zg["background"]
			B = np.sqrt(bkg["B_X"][0] ** 2 + bkg["B_Y"][0] ** 2
				+ bkg["B_Z"][0] ** 2)

			out = zg["output"]
			run["zones"][z] = {
				"nz": np.tensordot(w, out["nz"][:], axes=1),
				"counts": int(out["Nz"][:].sum()),
				"vol": vol,
				"Z_cen": Z_cen,
				"Z_corners": Zc,
				"B": B,
			}
	return run


def x_point_z(run, upper=False):
	"""
	Returns the height of the X-point (m), or None if there are no core
	zones.

	The core zones (those with a core face) cover the closed flux surfaces
	inside the separatrix, and the X-point is the separatrix's lowest point
	for a lower single null, or its highest for an upper single null. So
	it's taken as the lowest (or highest) corner of the core zones, which
	is within about a cell of the true X-point.

	Inputs
	 - run (dict)   : From load_run.
	 - upper (bool) : The divertor is at the top (upper single null).
	"""
	if not run["core_zones"]:
		return None
	if upper:
		return max(run["zones"][z]["Z_corners"].max()
			for z in run["core_zones"])
	return min(run["zones"][z]["Z_corners"].min() for z in run["core_zones"])


def print_summary(run, upper=False):
	"""
	Prints the source, redeposition, particle fates, mean time followed and
	inventory of a run.

	Inputs
	 - run (dict)    : From load_run.
	 - upper (bool)  : The divertor is at the top (upper single null).
	"""
	a = run["attrs"]
	S = a["source_rate"]
	print("Sputtered source:      {:.3e} atoms/s".format(S))
	if "redeposited_total" in a:
		redep = a["redeposited_total"]
		print("Redeposited:           {:.3e} atoms/s ({:.1f}% of the "
			"source)".format(redep, 100 * redep / S))
		print("  at wall faces:       {:.1f}%, into wall cells with no "
			"element: {:.1f}%".format(100 * a["redeposited_wall_face"] / S,
			100 * a["redeposited_unassigned"] / S))

	fates = {k[len("fates_"):]: int(v) for k, v in a.items()
		if k.startswith("fates_")}
	if fates:
		print("Particle fates (of {} particles):".format(run["imp_num"]))
		for k, v in fates.items():
			if v:
				print("  {:<12} {}".format(k, v))

	# Mean time followed: every recorded step is one time step of one
	# particle (assumes a constant time step)
	steps = sum(zd["counts"] for zd in run["zones"].values())
	t_mean = steps * run["dt"] / run["imp_num"]
	print("Mean time followed:    {:.3g} s ({:.3g} us)".format(t_mean,
		t_mean * 1e6))

	# Inventory, per zone and away from the divertor
	z_x = x_point_z(run, upper)
	total = 0.0
	away = 0.0
	per_zone = {}
	for z, zd in run["zones"].items():
		inv = zd["nz"] * zd["vol"]
		per_zone[z] = inv.sum()
		total += per_zone[z]
		if z_x is not None:
			above = zd["Z_cen"][:, :, None] > z_x
			away += inv[(above if not upper else ~above)
				& np.ones_like(inv, bool)].sum()
	print("Impurity inventory:    {:.3e} atoms, time-averaged (source x "
		"mean time followed = {:.3e})".format(total, S * t_mean))
	for z, inv in per_zone.items():
		if inv > 0:
			print("  zone {:>2}: {:.3e} ({:.1f}%)".format(z, inv,
				100 * inv / total if total else 0))
	if z_x is not None and total > 0:
		print("  {} the X-point (Z = {:.3f} m): {:.2f}%".format(
			"below" if upper else "above", z_x, 100 * away / total))


def wall_maps(run):
	"""
	Gross and net erosion per unit area at each place on the wall, summed
	over the toroidal cells.

	Inputs
	 - run (dict) : From load_run.

	Returns R and Z of each place (m), and the gross and net erosion per
	unit area there (atoms m^-2 s^-1). Elements at the same zone, cell
	and face in every toroidal cell are one place.
	"""
	key = np.stack([run["zone"], run["i"], run["j"], run["face"]], axis=1)
	_, first, inverse = np.unique(key, axis=0, return_index=True,
		return_inverse=True)
	inverse = inverse.ravel()
	area = np.bincount(inverse, weights=run["area"])
	gross = np.bincount(inverse, weights=run["source_mean"])
	net = np.bincount(inverse, weights=run.get("net_erosion",
		run["source_mean"]))
	with np.errstate(invalid="ignore", divide="ignore"):
		return (run["R_mid"][first], run["Z_mid"][first],
			np.where(area > 0, gross / area, 0.0),
			np.where(area > 0, net / area, 0.0))


def gyroperiod(run):
	"""
	Gyroperiod of a singly charged impurity ion at the plasma cell of the
	strongest source element (s), and that element's |B| (T).
	"""
	e = int(np.argmax(run["source_mean"]))
	zone, i, j, k = (int(run["zone"][e]), int(run["i"][e]),
		int(run["j"][e]), int(run["k"][e]))
	B = float(run["zones"][zone]["B"][i, j, k])
	m = run["mass_amu"] * AMU
	return 2 * np.pi * m / (EV * B), B


def plot_all(run, mesh_path=None, save=None):
	"""
	Plots the erosion maps, the time to landing and the charge at landing.

	Inputs
	 - run (dict)          : From load_run.
	 - mesh_path (str)     : SOLEDGE3X mesh file, for the wall outline.
	 - save (str or None)  : Prefix to save the figures to
	                         ({save}_erosion.png, {save}_landing.png). None
	                         shows them.
	"""
	walls = []
	if mesh_path is not None:
		import h5py
		with h5py.File(mesh_path, "r") as f:
			if "walls" in f:
				n = int(np.ravel(f["walls/Nwalls"][()])[0])
				for w in range(1, n + 1):
					g = f["walls/wall{}".format(w)]
					walls.append((g["R"][:], g["Z"][:]))

	# Gross and net erosion along the wall
	R, Z, gross, net = wall_maps(run)
	fig1, (ax1, ax2) = plt.subplots(1, 2, figsize=(11, 7))
	for ax in (ax1, ax2):
		for wR, wZ in walls:
			ax.plot(wR, wZ, color="k", lw=0.6)
		ax.scatter(R, Z, s=2, color="0.8", zorder=1)
		ax.set_aspect("equal")
		ax.set_xlabel("R (m)")
	ax1.set_ylabel("Z (m)")

	pos = gross > 0
	if pos.any():
		sc = ax1.scatter(R[pos], Z[pos], c=gross[pos], s=8, cmap="inferno",
			norm=LogNorm(vmin=gross[pos].max() * 1e-6, vmax=gross[pos].max()),
			zorder=2)
		fig1.colorbar(sc, ax=ax1, label="atoms m$^{-2}$ s$^{-1}$")
	ax1.set_title("Gross erosion")

	if "net_erosion" in run:
		lim = np.abs(net).max()
		if lim > 0:
			sc = ax2.scatter(R[net != 0], Z[net != 0], c=net[net != 0], s=8,
				cmap="RdBu_r", norm=SymLogNorm(linthresh=lim * 1e-4,
				vmin=-lim, vmax=lim), zorder=2)
			fig1.colorbar(sc, ax=ax2, label="atoms m$^{-2}$ s$^{-1}$")
		ax2.set_title("Net erosion (red) and deposition (blue)")
	fig1.tight_layout()

	# How quickly the redeposited atoms land, and their charge
	fig2, (ax3, ax4) = plt.subplots(1, 2, figsize=(11, 4.5))
	if "redeposited_lifetime" in run:
		edges = run["lifetime_edges"]
		hist = run["redeposited_lifetime"]
		total = hist.sum()
		if total > 0:
			# Fraction per decade, so the area under it is meaningful on a
			# log axis
			frac = hist / total / np.diff(np.log10(edges))
			ax3.stairs(frac, edges, fill=True, alpha=0.6)
			t_g, B = gyroperiod(run)
			ax3.axvline(t_g, color="k", ls="--",
				label="1+ gyroperiod at the\nstrongest source "
				"({:.2g} s, B = {:.2f} T)".format(t_g, B))
			ax3.legend(fontsize=8)
		ax3.set_xscale("log")
		ax3.set_xlabel("Time from being sputtered to landing (s)")
		ax3.set_ylabel("Fraction per decade")
		ax3.set_title("When the redeposited atoms land")

		q = run["redeposited_charge"]
		if q.sum() > 0:
			ax4.bar(np.arange(len(q)), q / q.sum())
			last = np.nonzero(q)[0].max()
			ax4.set_xlim(-0.5, last + 1.5)
		ax4.set_xlabel("Charge when landing")
		ax4.set_ylabel("Fraction")
		ax4.set_title("Charge of the redeposited atoms")
	fig2.tight_layout()

	if save is not None:
		fig1.savefig(save + "_erosion.png", dpi=150)
		fig2.savefig(save + "_landing.png", dpi=150)
		print("Saved {}_erosion.png and {}_landing.png".format(save, save))
	else:
		plt.show()


def main():
	parser = argparse.ArgumentParser(description="Summarize and plot a "
		"Flan run with imp_source_opt = \"sputter\".")
	parser.add_argument("nc_path", help="Flan output file")
	parser.add_argument("--mesh", default=None,
		help="SOLEDGE3X mesh file, to draw the wall outline")
	parser.add_argument("--upper", action="store_true",
		help="The divertor is at the top (upper single null)")
	parser.add_argument("--save", default=None,
		help="Save the figures with this prefix instead of showing them")
	parser.add_argument("--no-plots", action="store_true",
		help="Only print the summary")
	args = parser.parse_args()

	run = load_run(args.nc_path)
	print_summary(run, upper=args.upper)
	if not args.no_plots:
		plot_all(run, mesh_path=args.mesh, save=args.save)


if __name__ == "__main__":
	main()
