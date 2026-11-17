"""
Plots the impurity density (nz) from a Flan run on a SOLEDGE3X background
in the (R, Z) plane.

A SOLEDGE3X run's netCDF file has one group per zone (zone1, zone2, ...),
each with its own grid, so every zone is drawn on its own cell corners
(geometry/R_corners and geometry/Z_corners). The density can be shown for
one time frame or averaged over the frames, and for one toroidal plane or
averaged over them. Cells inside the wall (chi = 1) are shaded gray, and
the zone edges can be outlined. It can also make a movie of the time
frames, with one color scale for the whole movie.

Usage from the command line:

	python plot_s3x_nz.py s3x_example.nc
	python plot_s3x_nz.py s3x_example.nc --frame 0 --phi 3 --save nz.png
	python plot_s3x_nz.py s3x_example.nc --movie nz.mp4 --fps 10

or from Python:

	from plot_s3x_nz import plot_s3x_nz, movie_s3x_nz
	plot_s3x_nz("s3x_example.nc", frame=None, phi=None)
	movie_s3x_nz("s3x_example.nc", "nz.mp4", fps=10)

A movie ending in .mp4 needs ffmpeg; one ending in .gif is written with
Pillow, which comes with matplotlib.

This is a stopgap until the plotting tools handle multi-zone files.
"""

import argparse

import matplotlib.animation as animation
import matplotlib.pyplot as plt
from matplotlib.backends.backend_agg import FigureCanvasAgg
from matplotlib.figure import Figure
import netCDF4
import numpy as np
from matplotlib.colors import LogNorm, Normalize


def load_zones(nc, frame=None, phi=None):
	"""
	Loads each zone's cell corners, impurity density and wall mask.

	Inputs
	 - nc (netCDF4.Dataset) : Open Flan output file.
	 - frame (int, None     : Time frame to use, None to average over
	   or "all")              every frame, or "all" to keep every frame.
	 - phi (int or None)    : Toroidal plane (z index) to use, or None to
	                          average over every plane.

	Returns a list with a dictionary for each zone, holding the zone
	number ("zone"), the cell corners ("R", "Z", shape (nx+1, ny+1)), the
	density ("nz", shape (nx, ny), or (nt, nx, ny) with frame = "all") and
	the wall mask ("chi", shape (nx, ny), or None if the file doesn't have
	one).
	"""
	if "zones" not in nc.groups:
		raise ValueError("This file has a single zone. This script is only "
			"for multi-zone (e.g., SOLEDGE3X) runs.")

	nzones = len(nc["zones"]["nx"][:])
	zones = []
	for zone in range(1, nzones + 1):
		group = nc["zone{}".format(zone)]
		geo = group["geometry"]

		# Density is (time, x, y, z). Reduce over phi first (the last
		# axis), so only one plane per frame is kept in memory.
		var = group["output"]["nz"]
		if frame is None or frame == "all":
			nz = var[:].mean(axis=-1) if phi is None else var[:, :, :, phi]
			if frame is None:
				nz = nz.mean(axis=0)
		else:
			if not -var.shape[0] <= frame < var.shape[0]:
				raise IndexError("frame = {} is out of range: the file has "
					"{} time frames.".format(frame, var.shape[0]))
			nz = var[frame]
			nz = nz.mean(axis=-1) if phi is None else nz[:, :, phi]

		# The wall mask is axisymmetric, so any plane will do
		chi = geo["chi"][:, :, 0] if "chi" in geo.variables else None

		zones.append({"zone": zone, "R": geo["R_corners"][:],
			"Z": geo["Z_corners"][:], "nz": np.asarray(nz),
			"chi": None if chi is None else np.asarray(chi)})
	return zones


def make_norm(arrays, log=True, vmin=None, vmax=None):
	"""
	Color scale covering the densities in arrays (a list of numpy arrays).
	vmin and vmax default to the smallest nonzero (0 for a linear scale)
	and the largest density. Cells with no particles are left blank, so
	they don't set the scale.
	"""
	all_nz = np.concatenate([np.ravel(a) for a in arrays])
	positive = all_nz[all_nz > 0]
	if positive.size == 0:
		raise ValueError("Every cell has zero density.")
	if vmin is None:
		vmin = positive.min() if log else 0.0
	if vmax is None:
		vmax = positive.max()
	return LogNorm(vmin=vmin, vmax=vmax) if log \
		else Normalize(vmin=vmin, vmax=vmax)


def draw_zones(ax, zones, nz_list, norm, cmap="inferno", show_zones=True,
	show_wall=True):
	"""
	Draws the density of every zone on ax: the wall cells underneath (if
	show_wall), the density, and the zone edges (if show_zones).

	Inputs
	 - ax                   : Matplotlib axis.
	 - zones (list)         : From load_zones, for the corners and walls.
	 - nz_list (list)       : The density to draw in each zone, (nx, ny).
	 - norm                 : Color scale, from make_norm.

	Returns the density meshes, one per zone, so a movie can update them.
	"""
	meshes = []
	for z, nz in zip(zones, nz_list):
		R, Z = z["R"], z["Z"]

		# Wall cells underneath, so the density draws over them
		if show_wall and z["chi"] is not None:
			wall = np.ma.masked_where(z["chi"] < 0.5, z["chi"])
			ax.pcolormesh(R, Z, wall, cmap="Greys", vmin=0, vmax=2,
				shading="flat")

		meshes.append(ax.pcolormesh(R, Z, np.ma.masked_where(nz <= 0, nz),
			cmap=cmap, norm=norm, shading="flat"))

		if show_zones:
			edge_R = np.concatenate([R[0, :], R[1:, -1], R[-1, -2::-1],
				R[-2::-1, 0]])
			edge_Z = np.concatenate([Z[0, :], Z[1:, -1], Z[-1, -2::-1],
				Z[-2::-1, 0]])
			ax.plot(edge_R, edge_Z, color="k", lw=0.5)
	return meshes


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

	norm = make_norm([z["nz"] for z in zones], log, vmin, vmax)

	fig, ax = plt.subplots(figsize=(6, 8))
	mesh = draw_zones(ax, zones, [z["nz"] for z in zones], norm, cmap,
		show_zones, show_wall)[-1]

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


def movie_s3x_nz(nc_path, out, phi=None, start=0, stop=None, step=1,
	fps=10, log=True, vmin=None, vmax=None, cmap="inferno", show_zones=True,
	show_wall=True, dpi=150):
	"""
	Makes a movie of the impurity density of a multi-zone (SOLEDGE3X) Flan
	run in the (R, Z) plane, one movie frame per time frame.

	The color scale is the same for the whole movie, so changes in density
	from frame to frame show. By default it covers the smallest nonzero to
	the largest density over the frames in the movie; set vmin and vmax to
	choose it.

	Inputs
	 - nc_path (str)        : Path to the Flan netCDF file.
	 - out (str)            : Movie file to write. .mp4 needs ffmpeg; .gif
	                          is written with Pillow.
	 - phi (int or None)    : Toroidal plane (z index) to show, or None for
	                          the toroidal average.
	 - start, stop, step    : Which time frames to use, as in
	   (int)                  range(start, stop, step). stop = None goes to
	                          the last frame.
	 - fps (float)          : Movie frames per second.
	 - log (bool)           : Log color scale if True, linear if False.
	 - vmin, vmax (float)   : Color scale limits.
	 - cmap (str)           : Matplotlib colormap name.
	 - show_zones (bool)    : Outline the edges of each zone.
	 - show_wall (bool)     : Shade cells inside the wall (chi = 1) gray.
	 - dpi (int)            : Resolution of the movie frames.

	Returns the path of the movie written.
	"""
	with netCDF4.Dataset(nc_path) as nc:
		zones = load_zones(nc, "all", phi)
		times = np.asarray(nc["geometry"]["time"][:])

	frames = list(range(len(times)))[start:stop:step]
	if not frames:
		raise ValueError("No time frames selected (the file has {}, "
			"start = {}, stop = {}, step = {}).".format(len(times), start,
			stop, step))

	# One color scale over every frame in the movie
	norm = make_norm([z["nz"][frames] for z in zones], log, vmin, vmax)

	# Drawn off screen (no pyplot window), so it works without a display
	# and doesn't change the backend of an interactive session
	fig = Figure(figsize=(6, 8), layout="constrained")
	FigureCanvasAgg(fig)
	ax = fig.add_subplot()
	meshes = draw_zones(ax, zones, [z["nz"][frames[0]] for z in zones],
		norm, cmap, show_zones, show_wall)
	cbar = fig.colorbar(meshes[-1], ax=ax)
	cbar.set_label(r"$n_z$ (m$^{-3}$)", fontsize=14)
	where = "toroidal average" if phi is None \
		else "toroidal plane {}".format(phi)
	title = ax.set_title("", fontsize=12)
	ax.set_aspect("equal")
	ax.set_xlabel("R (m)", fontsize=14)
	ax.set_ylabel("Z (m)", fontsize=14)

	def update(i):
		f = frames[i]
		for mesh, z in zip(meshes, zones):
			nz = z["nz"][f]
			mesh.set_array(np.ma.masked_where(nz <= 0, nz))
		title.set_text("t = {:.6g} s (frame {}), {}".format(times[f], f,
			where))
		return meshes + [title]

	anim = animation.FuncAnimation(fig, update, frames=len(frames),
		blit=False)
	if str(out).lower().endswith(".gif"):
		writer = animation.PillowWriter(fps=fps)
	else:
		if not animation.writers.is_available("ffmpeg"):
			raise RuntimeError("Writing {} needs ffmpeg, which wasn't "
				"found. Install it, or write a .gif instead.".format(out))
		writer = animation.FFMpegWriter(fps=fps)
	anim.save(out, writer=writer, dpi=dpi)
	print("Wrote {} ({} frames)".format(out, len(frames)))
	return out


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
	parser.add_argument("--movie", default=None,
		help="make a movie of the time frames and write it to this path "
		"(.mp4 needs ffmpeg, .gif doesn't); --frame and --save are then "
		"ignored")
	parser.add_argument("--fps", type=float, default=10,
		help="movie frames per second (default 10)")
	parser.add_argument("--start", type=int, default=0,
		help="first time frame of the movie (default 0)")
	parser.add_argument("--stop", type=int, default=None,
		help="time frame the movie stops before (default: to the end)")
	parser.add_argument("--step", type=int, default=1,
		help="use every step-th time frame in the movie (default 1)")
	parser.add_argument("--dpi", type=int, default=150,
		help="movie resolution (default 150)")
	args = parser.parse_args()

	if args.movie is not None:
		movie_s3x_nz(args.nc_path, args.movie, phi=args.phi,
			start=args.start, stop=args.stop, step=args.step, fps=args.fps,
			log=not args.linear, vmin=args.vmin, vmax=args.vmax,
			cmap=args.cmap, show_zones=not args.no_zones,
			show_wall=not args.no_wall, dpi=args.dpi)
		return

	# Nothing is shown, so draw off screen: this works without a display
	if args.no_show:
		plt.switch_backend("Agg")

	plot_s3x_nz(args.nc_path, frame=args.frame, phi=args.phi,
		log=not args.linear, vmin=args.vmin, vmax=args.vmax, cmap=args.cmap,
		show_zones=not args.no_zones, show_wall=not args.no_wall,
		save=args.save, showplot=not args.no_show)


if __name__ == "__main__":
	main()
