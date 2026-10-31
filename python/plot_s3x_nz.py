"""
Plots the impurity density (nz) from a Flan run on a SOLEDGE3X background
in the (R, Z) plane.

A SOLEDGE3X run's netCDF file has one group per zone (zone1, zone2, ...),
each with its own grid, so every zone is drawn on its own cell corners
(geometry/R_corners and geometry/Z_corners). The density can be shown for
one time frame or averaged over the frames, and for one toroidal plane or
averaged over them. Cells inside the wall (chi = 1) are shaded gray, and
the zone edges can be outlined.

Usage from the command line:

	python plot_s3x_nz.py s3x_example.nc
	python plot_s3x_nz.py s3x_example.nc --frame 0 --phi 3 --save nz.png

or from Python:

	from plot_s3x_nz import plot_s3x_nz
	plot_s3x_nz("s3x_example.nc", frame=None, phi=None)

This is a stopgap until the plotting tools handle multi-zone files.
"""

import argparse

import matplotlib.pyplot as plt
import netCDF4
import numpy as np
from matplotlib.colors import LogNorm, Normalize


def load_zones(nc, frame=None, phi=None):
	"""
	Loads each zone's cell corners, impurity density and wall mask.

	Inputs
	 - nc (netCDF4.Dataset) : Open Flan output file.
	 - frame (int or None)  : Time frame to use, or None to average over
	                          every frame.
	 - phi (int or None)    : Toroidal plane (z index) to use, or None to
	                          average over every plane.

	Returns a list with a dictionary for each zone, holding the zone
	number ("zone"), the cell corners ("R", "Z", shape (nx+1, ny+1)), the
	density ("nz", shape (nx, ny)) and the wall mask ("chi", shape (nx, ny),
	or None if the file doesn't have one).
	"""
	if "zones" not in nc.groups:
		raise ValueError("This file has a single zone. This script is only "
			"for multi-zone (e.g., SOLEDGE3X) runs.")

	nzones = len(nc["zones"]["nx"][:])
	zones = []
	for zone in range(1, nzones + 1):
		group = nc["zone{}".format(zone)]
		geo = group["geometry"]

		# Density is (time, x, y, z)
		nz = group["output"]["nz"][:]
		nz = nz.mean(axis=0) if frame is None else nz[frame]
		nz = nz.mean(axis=2) if phi is None else nz[:, :, phi]

		# The wall mask is axisymmetric, so any plane will do
		chi = geo["chi"][:, :, 0] if "chi" in geo.variables else None

		zones.append({"zone": zone, "R": geo["R_corners"][:],
			"Z": geo["Z_corners"][:], "nz": np.asarray(nz),
			"chi": None if chi is None else np.asarray(chi)})
	return zones


def plot_s3x_nz(nc_path, frame=None, phi=None, log=True, vmin=None,
	vmax=None, cmap="inferno", show_zones=True, show_wall=True, save=None,
	showplot=True):
	"""
	Plots the impurity density of a multi-zone (SOLEDGE3X) Flan run in the
	(R, Z) plane.

	Inputs
	 - nc_path (str)        : Path to the Flan netCDF file.
	 - frame (int or None)  : Time frame to plot, or None to average over
	                          every frame.
	 - phi (int or None)    : Toroidal plane (z index) to plot, or None to
	                          average over every plane.
	 - log (bool)           : Log color scale if True, linear if False.
	 - vmin, vmax (float)   : Color scale limits. Default to the smallest
	                          nonzero and largest density.
	 - cmap (str)           : Matplotlib colormap name.
	 - show_zones (bool)    : Outline the edges of each zone.
	 - show_wall (bool)     : Shade cells inside the wall (chi = 1) gray.
	 - save (str or None)   : Path to save the figure to, if given.
	 - showplot (bool)      : Show the figure.

	Returns the matplotlib figure and axis.
	"""
	with netCDF4.Dataset(nc_path) as nc:
		zones = load_zones(nc, frame, phi)
		times = nc["geometry"]["time"][:]

	# Color scale from every zone. Cells with no particles are left blank.
	all_nz = np.concatenate([z["nz"].ravel() for z in zones])
	positive = all_nz[all_nz > 0]
	if positive.size == 0:
		raise ValueError("Every cell has zero density.")
	if vmin is None:
		vmin = positive.min() if log else 0.0
	if vmax is None:
		vmax = positive.max()
	norm = LogNorm(vmin=vmin, vmax=vmax) if log \
		else Normalize(vmin=vmin, vmax=vmax)

	fig, ax = plt.subplots(figsize=(6, 8))
	mesh = None
	for z in zones:
		R, Z = z["R"], z["Z"]

		# Wall cells underneath, so the density draws over them
		if show_wall and z["chi"] is not None:
			wall = np.ma.masked_where(z["chi"] < 0.5, z["chi"])
			ax.pcolormesh(R, Z, wall, cmap="Greys", vmin=0, vmax=2,
				shading="flat")

		nz = np.ma.masked_where(z["nz"] <= 0, z["nz"])
		mesh = ax.pcolormesh(R, Z, nz, cmap=cmap, norm=norm,
			shading="flat")

		if show_zones:
			edge_R = np.concatenate([R[0, :], R[1:, -1], R[-1, -2::-1],
				R[-2::-1, 0]])
			edge_Z = np.concatenate([Z[0, :], Z[1:, -1], Z[-1, -2::-1],
				Z[-2::-1, 0]])
			ax.plot(edge_R, edge_Z, color="k", lw=0.5)

	cbar = fig.colorbar(mesh, ax=ax)
	cbar.set_label(r"$n_z$ (m$^{-3}$)", fontsize=14)

	# What's being shown
	when = "average of {} frames".format(len(times)) if frame is None \
		else "t = {:.6g} s".format(times[frame])
	where = "toroidal average" if phi is None \
		else "toroidal plane {}".format(phi)
	ax.set_title("{}, {}".format(when, where), fontsize=12)

	ax.set_aspect("equal")
	ax.set_xlabel("R (m)", fontsize=14)
	ax.set_ylabel("Z (m)", fontsize=14)
	fig.tight_layout()

	if save is not None:
		fig.savefig(save, dpi=200)
	if showplot:
		plt.show()

	return fig, ax


def main():
	parser = argparse.ArgumentParser(description="Plot the impurity "
		"density of a Flan run on a SOLEDGE3X background in the (R, Z) "
		"plane.")
	parser.add_argument("nc_path", help="Flan netCDF file")
	parser.add_argument("--frame", type=int, default=None,
		help="time frame to plot (default: average over frames)")
	parser.add_argument("--phi", type=int, default=None,
		help="toroidal plane to plot (default: toroidal average)")
	parser.add_argument("--linear", action="store_true",
		help="linear color scale instead of log")
	parser.add_argument("--vmin", type=float, default=None,
		help="color scale minimum")
	parser.add_argument("--vmax", type=float, default=None,
		help="color scale maximum")
	parser.add_argument("--cmap", default="inferno", help="colormap")
	parser.add_argument("--no-zones", action="store_true",
		help="don't outline the zones")
	parser.add_argument("--no-wall", action="store_true",
		help="don't shade the wall cells")
	parser.add_argument("--save", default=None,
		help="save the figure to this path")
	parser.add_argument("--no-show", action="store_true",
		help="don't open a window (e.g., when only saving)")
	args = parser.parse_args()

	plot_s3x_nz(args.nc_path, frame=args.frame, phi=args.phi,
		log=not args.linear, vmin=args.vmin, vmax=args.vmax, cmap=args.cmap,
		show_zones=not args.no_zones, show_wall=not args.no_wall,
		save=args.save, showplot=not args.no_show)


if __name__ == "__main__":
	main()
