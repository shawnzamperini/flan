"""
Checks the conventions of a JAX-Hermes3 run before it is used as a Flan
background.

A JAX-Hermes3 run is a hypnotoad grid (netCDF) plus one .npz snapshot per
time with Ne, Pe, Pi, NVi, NVe, phi and Vort on the (x, y, z) cell centers,
all normalized (Nnorm, Tnorm, Bnorm). The case's JSON file gives the
normalization, the species masses, the conduction coefficients and
zperiod. Flan's HERMES reader relies on several conventions that neither
file states, so this script tests each one on the data:

1. Cell corners. hypnotoad stores the four corners of each cell. The
   script finds which corner is which by matching neighboring cells, and
   checks that cells on either side of the branch cuts (core and private
   flux region) share their corners.
2. Handedness. The grid's Jacobian J and curvature bxcvx are compared with
   those computed from Rxy, Zxy and Btxy, assuming (R, z, Z) is right
   handed (z, the toroidal angle, increasing counterclockwise seen from
   above). From that the direction of the ion grad-B drift is found and
   compared with the X-point position.
3. Field alignment. Between neighboring y planes, the toroidal shift at
   which the density correlates best is compared with the change in
   zShift (plus ShiftAngle across the core branch cut). This sets the
   sign used for the field-aligned parallel derivatives below.
4. Parallel flow. V = NV / (AA * N) should be at least sonic and point
   into the targets: negative at the first y cell, positive at the last.
5. Particle and energy balance. The particle and energy fluxes through the
   y faces next to both targets are summed and compared with the sources
   in the case's input (ne_flux, and heating_W per species). The parallel
   heat flux is q = -kappa grad_par T, with kappa as in Hermes-3's
   braginskii_conduction (kappa = c P tau / m, tau = 1 / sum of the
   species' Braginskii collision frequencies), and the gradient taken
   along the field line. The same balance with the gradient at fixed z
   shows how wrong a plain y-derivative is on this grid.

Usage:

	python check_hermes_conventions.py --root /path/to/data \\
		--case my_case --grid mesh/reversed.nc

The forward-field grid is used for the grad-B drift test as well if
--grid-other is given (e.g. mesh/forward.nc).
"""

import argparse
import glob
import json
import os
import re

import netCDF4
import numpy as np

# SI constants, as in BOUT++ (SI::qe, SI::Me, SI::Mp, SI::e0).
QE = 1.602176634e-19
ME = 9.1093837015e-31
MP = 1.67262192369e-27
E0 = 8.8541878128e-12

# Number of x guard cells on each side of a BOUT++ grid.
MXG = 2


def soft_floor(value, minval):
	"""
	Hermes-3's softFloor: values below zero are set to zero, then a smooth
	floor is added so the result is never below minval.

	@param value Array to floor.
	@param minval Floor (> 0).
	@return The floored array.
	"""
	v = np.maximum(value, 0.0)
	return v + minval * np.exp(-v / minval)


class Grid:
	"""
	The parts of a hypnotoad grid file used here, with the single-null
	topology (where the branch cuts are) worked out from the separatrix
	indices.
	"""

	def __init__(self, path):
		"""
		@param path Path to the hypnotoad netCDF grid file.
		"""
		with netCDF4.Dataset(path) as d:
			self.v = {k: np.array(d[k][:]) for k in d.variables
				if d[k].dtype != str and d[k].ndim <= 2}
		self.nx = int(self.v["nx"])
		self.ny = int(self.v["ny"])
		self.ixseps1 = int(self.v["ixseps1"])
		self.j11 = int(self.v["jyseps1_1"])
		self.j22 = int(self.v["jyseps2_2"])
		if int(self.v["jyseps2_1"]) != int(self.v["jyseps1_2"]):
			raise SystemExit(f"{path} is not a single null "
				"(jyseps2_1 != jyseps1_2); only single null is handled.")
		if int(self.v["y_boundary_guards"]) != 0:
			raise SystemExit(f"{path} has y guard cells; not handled.")

	def __getitem__(self, key):
		return self.v[key]

	def y_next(self, i, j):
		"""
		The cell after (i, j) in +y along the field, following the branch
		cuts.

		For x < ixseps1 (core and private flux region) the core's last
		cell (jyseps2_2) connects to its first (jyseps1_1 + 1), and the
		inner leg's last cell (jyseps1_1) to the outer leg's first
		(jyseps2_2 + 1). In the SOL the y index runs from target to target.

		@param i x index.
		@param j y index.
		@return (jn, cut): the next y index (None at the outer target) and
		        whether the step crosses the core branch cut, where
		        ShiftAngle has to be added to the change in zShift.
		"""
		if i < self.ixseps1:
			if j == self.j22:
				return self.j11 + 1, True
			if j == self.j11:
				return self.j22 + 1, False
		if j == self.ny - 1:
			return None, False
		return j + 1, False


def load_case(root, case):
	"""
	Reads a case's JSON file and lists its snapshots in time order.

	@param root Root of the JAX-Hermes3 data (holds inputs/ and data/).
	@param case Case name (the JSON file is inputs/<case>.json).
	@return (info, snapshot paths).
	"""
	with open(os.path.join(root, "inputs", case + ".json")) as f:
		info = json.load(f)
	snaps = sorted(glob.glob(os.path.join(root, "data", info["resolution"],
		case, "snapshots", "snap*.npz")))
	if not snaps:
		raise SystemExit(f"No snapshots for {case} under {root}/data.")
	return info, snaps


def parse_norm(info):
	"""
	Pulls the normalization out of the case's "norm" string, which looks
	like "Norm(Nnorm=1e+19, Tnorm=50.0, Bnorm=1.0)".

	@param info The case's JSON contents.
	@return Dict with Nnorm, Tnorm and Bnorm.
	"""
	return {k: float(v) for k, v in
		re.findall(r"(\w+)=([-+.\deE]+)", info["pars"]["norm"])}


def parse_species(text):
	"""
	Reads AA and charge from a "Species(AA=2.0, charge=1.0)" string.

	@param text The string from the case's pars.
	@return (AA, charge).
	"""
	d = {k: float(v) for k, v in re.findall(r"(\w+)=([-+.\deE]+)", text)}
	return d["AA"], d["charge"]


# ---------------------------------------------------------------------------
# 1. Corners
# ---------------------------------------------------------------------------

def check_corners(g):
	"""
	Works out which stored corner is which and checks that neighbors and
	cells across the branch cuts share corners.

	hypnotoad stores Rxy_corners plus the lower_right, upper_right and
	upper_left corners. Whether "right" means +x or +y isn't stated, so
	both readings are tested on neighboring cells and the one that matches
	is reported.

	@param g The grid.
	"""
	print("\n== 1. Cell corners ==")
	names = {"ll": "corners", "lr": "lower_right_corners",
		"ur": "upper_right_corners", "ul": "upper_left_corners"}
	R = {k: g["Rxy_" + v] for k, v in names.items()}
	Z = {k: g["Zxy_" + v] for k, v in names.items()}
	size = np.hypot(R["ur"] - R["ll"], Z["ur"] - Z["ll"])

	def mismatch(a, ia, b, ib):
		# Distance between corner a of cells ia and corner b of cells ib,
		# relative to the cell size.
		return np.hypot(R[a][ia] - R[b][ib], Z[a][ia] - Z[b][ib]) / size[ia]

	# Neighbors in +x (all y) and +y (only within one region, no cut).
	ix = (slice(0, g.nx - 1), slice(None))
	ixp = (slice(1, g.nx), slice(None))
	jlist = [j for j in range(g.ny - 1) if j not in (g.j11, g.j22)]
	iy = (slice(None), jlist)
	iyp = (slice(None), [j + 1 for j in jlist])

	readings = {
		"right = +x, upper = +y": (("lr", "ll"), ("ul", "ll")),
		"right = +y, upper = +x": (("ul", "ll"), ("lr", "ll")),
	}
	best = None
	for label, ((ax, bx), (ay, by)) in readings.items():
		mx = max(mismatch(ax, ix, bx, ixp).max(),
			mismatch(ay, iy, by, iyp).max())
		print(f"  {label}: max neighbor mismatch {mx:.2e} cell sizes")
		if best is None or mx < best[1]:
			best = (label, mx, (ax, ay))
	label, mx, (ax, ay) = best
	print(f"  -> corners follow '{label}'")

	# Across the cuts, for x < ixseps1: the core's last cell's +y corners
	# against the core's first cell's -y corners, same for the legs.
	core = np.arange(0, g.ixseps1)
	for name, ja, jb in (("core cut", g.j22, g.j11 + 1),
			("PFR cut", g.j11, g.j22 + 1)):
		m = mismatch(ay, (core, ja), "ll", (core, jb)).max()
		print(f"  {name} (y {ja} -> {jb}): max mismatch {m:.2e} "
			"cell sizes")


# ---------------------------------------------------------------------------
# 2. Handedness
# ---------------------------------------------------------------------------

def basis(g):
	"""
	Covariant basis vectors e_x = dr/dx and e_y = dr/dy in the poloidal
	plane, from centered differences of Rxy and Zxy in the index space,
	divided by dx and dy. Only interior cells away from the branch cuts are
	filled; the rest are NaN.

	@param g The grid.
	@return (aR, aZ, bR, bZ): R and Z parts of e_x and e_y.
	"""
	R, Z = g["Rxy"], g["Zxy"]
	shape = R.shape
	aR, aZ, bR, bZ = (np.full(shape, np.nan) for _ in range(4))
	aR[1:-1] = (R[2:] - R[:-2]) / (2 * g["dx"][1:-1])
	aZ[1:-1] = (Z[2:] - Z[:-2]) / (2 * g["dx"][1:-1])
	for j in range(1, g.ny - 1):
		if j in (g.j11, g.j11 + 1, g.j22, g.j22 + 1):
			continue
		bR[:, j] = (R[:, j + 1] - R[:, j - 1]) / (2 * g["dy"][:, j])
		bZ[:, j] = (Z[:, j + 1] - Z[:, j - 1]) / (2 * g["dy"][:, j])
	# The x difference is invalid next to the separatrix only where the
	# regions on either side differ, which isn't the case for a single null
	# along x, so nothing more is masked.
	return aR, aZ, bR, bZ


def check_handedness(g, label):
	"""
	Tests whether the grid's toroidal angle z makes (R, z, Z) right handed,
	in two independent ways, then reports the ion grad-B drift direction.

	For a right-handed (R, phi, Z) with z = phi, e_z = R phi_hat and

	    J = e_x . (e_y x e_z) = R (aZ bR - aR bZ),

	where (aR, aZ) = e_x and (bR, bZ) = e_y. Components of e_x and e_y
	along phi_hat (from the shifted metric) don't change the product.
	For a left-handed one the sign flips. The curvature test uses
	b x kappa ~ (Bt/B) phi_hat x (-R_hat/R) = (Bt / (B R)) Z_hat (right
	handed), whose x component has the sign of Bt * (grad x)_Z.

	@param g The grid.
	@param label Name of the grid for the printout.
	@return +1 if right handed, -1 if left handed.
	"""
	print(f"\n== 2. Handedness: {label} ==")
	aR, aZ, bR, bZ = basis(g)
	R = g["Rxy"]
	J_rh = R * (aZ * bR - aR * bZ)
	ok = np.isfinite(J_rh)
	ok[:MXG] = ok[-MXG:] = False
	ratio = g["J"][ok] / J_rh[ok]
	agree_J = np.mean(ratio > 0)
	print(f"  J(grid) / J(right handed): median {np.median(ratio):+.3f}, "
		f"same sign in {100 * agree_J:.1f}% of {ok.sum()} cells")

	# Curvature: (grad x)_Z from inverting [[aR, bR], [aZ, bZ]].
	det = aR * bZ - aZ * bR
	gx_R, gx_Z = bZ / det, -bR / det
	steep = ok & (np.abs(gx_Z) > 0.3 * np.hypot(gx_R, gx_Z))
	pred = np.sign(g["Btxy"] * gx_Z)
	agree_c = np.mean(np.sign(g["bxcvx"][steep]) == pred[steep])
	print(f"  sign(bxcvx) matches right-handed b x kappa in "
		f"{100 * agree_c:.1f}% of {steep.sum()} cells")

	hand = 1 if agree_J > 0.5 else -1
	if (agree_c > 0.5) != (agree_J > 0.5):
		print("  !! The two tests disagree.")
	print(f"  -> (R, z, Z) is {'right' if hand > 0 else 'left'} handed")

	# Ion grad-B drift ~ B x grad B. At the outboard midplane grad B points
	# in -R, so (B x grad B)_Z = -B_phi dB/dR, with B_phi = hand * Btxy.
	isep = g.ixseps1
	jmid = g.j11 + 1 + np.argmax(R[isep, g.j11 + 1:g.j22 + 1])
	B = g["Bxy"][:, jmid]
	dBdR = (B[isep + 1] - B[isep - 1]) / (R[isep + 1, jmid]
		- R[isep - 1, jmid])
	drift_Z = -hand * g["Btxy"][isep, jmid] * dBdR
	zx = 0.25 * (g["Zxy_upper_right_corners"][isep - 1, g.j11]
		+ g["Zxy_upper_left_corners"][isep, g.j11]
		+ g["Zxy_upper_right_corners"][isep - 1, g.j22]
		+ g["Zxy_upper_left_corners"][isep, g.j22])
	zaxis = np.mean(g["Zxy"][MXG, g.j11 + 1:g.j22 + 1])
	toward = np.sign(drift_Z) == np.sign(zx - zaxis)
	print(f"  Btxy at the outboard midplane {g['Btxy'][isep, jmid]:+.3f} T; "
		f"X-point near Z = {zx:+.3f} m (axis Z ~ {zaxis:+.3f} m)")
	print(f"  -> ion grad-B drift points {'down' if drift_Z < 0 else 'up'}, "
		f"{'toward' if toward else 'away from'} the X-point")
	return hand


# ---------------------------------------------------------------------------
# 3. Field alignment
# ---------------------------------------------------------------------------

def y_shift_cells(g, i, j, jn, cut, dz):
	"""
	Toroidal shift, in z cells, of a field line going from plane j to jn.

	@param g The grid.
	@param i x index.
	@param j, jn y indices of the two planes.
	@param cut Whether the step crosses the core branch cut.
	@param dz Toroidal cell width (rad).
	@return zShift(i, jn) - zShift(i, j) (+ ShiftAngle across the cut),
	        in cells.
	"""
	d = g["zShift"][i, jn] - g["zShift"][i, j]
	if cut:
		d += g["ShiftAngle"][i]
	return d / dz


def shift_z(f, s):
	"""
	Samples a periodic array at k + s along its last axis, with linear
	interpolation.

	@param f Array (..., nz).
	@param s Shift in cells (scalar, or broadcastable to f[..., 0]).
	@return f at k + s for each k.
	"""
	nz = f.shape[-1]
	k = np.arange(nz) + np.asarray(s)[..., None]
	k0 = np.floor(k).astype(int)
	w = k - k0
	return ((1 - w) * np.take_along_axis(f, k0 % nz, axis=-1)
		+ w * np.take_along_axis(f, (k0 + 1) % nz, axis=-1))


def check_alignment(g, ne, dz):
	"""
	Finds, for a sample of neighboring y planes, the toroidal shift at
	which the density correlates best, and compares it with +/- the
	change in zShift.

	@param g The grid.
	@param ne Density (nx, ny, nz) from one snapshot.
	@param dz Toroidal cell width (rad).
	@return +1 if f(jn) at z + shift lines up with f(j) at z, else -1.
	"""
	print("\n== 3. Field alignment ==")
	nz = ne.shape[2]
	votes = {+1: 0, -1: 0}
	rows = []
	isamp = [g.ixseps1 - 6, g.ixseps1 + 3, g.ixseps1 + 15]
	for i in isamp:
		for j in range(g.ny):
			jn, cut = g.y_next(i, j)
			if jn is None or (j % 4 and not cut):
				continue
			a = ne[i, j] - ne[i, j].mean()
			b = ne[i, jn] - ne[i, jn].mean()
			# corr[m] = sum_k a[k] b[k + m]
			corr = np.fft.ifft(np.conj(np.fft.fft(a)) * np.fft.fft(b)).real
			m = int(np.argmax(corr))
			s = y_shift_cells(g, i, j, jn, cut, dz)
			d = {sg: abs((m - sg * s + nz / 2) % nz - nz / 2)
				for sg in (+1, -1)}
			sg = min(d, key=d.get)
			votes[sg] += 1
			rows.append((i, j, jn, cut, s, m, d[sg]))
	for i, j, jn, cut, s, m, err in rows:
		print(f"  x {i:2d}, y {j:2d} -> {jn:2d}{' (cut)' if cut else '      '}: "
			f"shift {s:+8.1f} cells ({s % nz:5.1f} mod nz), best "
			f"correlation at {m:2d}, off by {err:.1f}")
	sign = max(votes, key=votes.get)
	print(f"  -> f(y+1) at z {'+' if sign > 0 else '-'} shift matches f(y) "
		f"at z in {votes[sign]} of {sum(votes.values())} pairs")
	return sign


# ---------------------------------------------------------------------------
# 4.-5. Flows, collisions, conduction and balances
# ---------------------------------------------------------------------------

def collision_frequencies(ne, te, ti, ai):
	"""
	Total Braginskii collision frequencies of the electrons and the ions,
	as in Hermes-3's braginskii_collisions with its defaults: e-e, e-i and
	i-i collisions, with i-e from momentum conservation,
	m_e n_e nu_ei = m_i n_i nu_ie, and a single ion species of charge 1
	with n_i = n_e.

	@param ne Density (m^-3).
	@param te, ti Temperatures (eV).
	@param ai Ion mass number.
	@return (nu_e, nu_i) in s^-1.
	"""
	# e-e (NRL formulary 2019, p. 34).
	telim = soft_floor(te, 0.1)
	nelim = soft_floor(ne, 1e10)
	logte = np.log(telim)
	clog = 30.4 - 0.5 * np.log(nelim) + 1.25 * logte \
		- np.sqrt(1e-5 + (logte - 2) ** 2 / 16)
	v1sq = 2 * telim * QE / ME
	nu_ee = QE ** 4 * np.maximum(ne, 0) * soft_floor(clog, 1.0) * 2 \
		/ (3 * (np.pi * 2 * v1sq) ** 1.5 * (E0 * ME) ** 2)

	# e-i, with the Coulomb log branches of Hermes-3. Zi = 1, ni = ne.
	me_mi = ME / (MP * ai)
	with np.errstate(divide="ignore", invalid="ignore"):
		clog = np.where((te < 0.1) | (ne < 1e10), 10.0,
			np.where(te < ti * me_mi,
				23 - 0.5 * np.log(ne) + 1.5 * np.log(ti) - np.log(ai),
				np.where(te < np.exp(2),
					30.0 - 0.5 * np.log(ne) + 1.5 * np.log(te),
					31.0 - 0.5 * np.log(ne) + np.log(te))))
	vesq = 2 * soft_floor(te, 0.1) * QE / ME
	visq = 2 * soft_floor(ti, 0.1) * QE / (MP * ai)
	nu_ei = QE ** 4 * np.maximum(ne, 0) * soft_floor(clog, 1.0) \
		* (1 + me_mi) / (3 * (np.pi * (vesq + visq)) ** 1.5 * (E0 * ME) ** 2)

	# i-i, same species.
	tlim = soft_floor(ti, 0.1)
	nlim = soft_floor(ne, 1e10)
	clog = 29.91 - np.log((2 * ai) / (2 * ai * tlim) * np.sqrt(2 * nlim / tlim))
	vsq = 2 * tlim * QE / (MP * ai)
	nu_ii = QE ** 4 * nlim * soft_floor(clog, 1.0) * 2 \
		/ (3 * (np.pi * 2 * vsq) ** 1.5 * (E0 * MP * ai) ** 2)

	# i-e from m_e n_e nu_ei = m_i n_i nu_ie (Hermes-3 uses the electron
	# AA = 1/1836 here, so m_e/m_i = (1/1836)/ai).
	nu_ie = nu_ei * (1 / 1836) / ai
	return nu_ee + nu_ei, nu_ii + nu_ie


def load_snapshot(path, norm, ai, ae):
	"""
	Reads one snapshot and converts it to SI units.

	@param path Path to the .npz file.
	@param norm Dict with Nnorm and Tnorm.
	@param ai, ae Ion and electron mass numbers.
	@return Dict with ne (m^-3), te, ti (eV), vi, ve (m/s), all (nx, ny, nz).
	"""
	z = np.load(path)
	cs0 = np.sqrt(QE * norm["Tnorm"] / MP)
	n = z["Ne"]
	nf = soft_floor(n, 1e-5)
	return {
		"t": float(z["t"]),
		"ne": n * norm["Nnorm"],
		"te": z["Pe"] / nf * norm["Tnorm"],
		"ti": z["Pi"] / nf * norm["Tnorm"],
		"vi": z["NVi"] / (ai * nf) * cs0,
		"ve": z["NVe"] / (ae * nf) * cs0,
	}


def check_flow(g, s, ai):
	"""
	Prints the parallel ion flow at the first and last y cells of the SOL,
	in units of the sound speed sqrt((Te + Ti) / m_i).

	@param g The grid.
	@param s A snapshot in SI units (load_snapshot).
	@param ai Ion mass number.
	"""
	print("\n== 4. Parallel flow at the targets ==")
	isol = slice(g.ixseps1, g.nx - MXG)
	for name, j in (("first y cell (inner target)", 0),
			("last y cell (outer target)", g.ny - 1)):
		cs = np.sqrt((s["te"][isol, j] + s["ti"][isol, j]) * QE / (ai * MP))
		m = s["vi"][isol, j] / cs
		print(f"  {name}: V/cs median {np.median(m):+.2f}, "
			f"5-95% [{np.percentile(m, 5):+.2f}, {np.percentile(m, 95):+.2f}]")


def face_fluxes(g, s, info, j, align, dz):
	"""
	Particle flux and parallel energy flux (W) through the y face between
	cells j and j + 1, summed over x (physical cells) and the whole torus,
	positive toward +y.

	The face area normal to b is J dx dz / sqrt(g_22) at the face (the
	_ylow fields of cell j + 1), times zperiod for the full torus. Face
	values are the average of the two cells; the temperature gradient is
	(T(j + 1) - T(j)) / L, with L = (sqrt(g_22) dy) averaged over the two
	cells. With align = +/-1, the j + 1 values are taken along the field
	line (z +/- the zShift change); with align = 0 at fixed z.

	The energy flux per species is 5/2 n T V + q, plus 1/2 m_i n V^3 for
	the ions, with q = -kappa grad_par T and kappa from
	braginskii_conduction.

	@param g The grid.
	@param s A snapshot in SI units.
	@param info The case's JSON contents.
	@param j y index of the cell before the face.
	@param align +1, -1 or 0, as above.
	@param dz Toroidal cell width (rad).
	@return Dict with particles (s^-1), and conducted and total energy
	        fluxes (W) for each species.
	"""
	pars = info["pars"]
	ai, _ = parse_species(pars["ion"])
	ae, _ = parse_species(pars["electron"])
	x = slice(MXG, g.nx - MXG)
	jn = j + 1

	def pair(f):
		a = f[x, j]
		b = f[x, jn]
		if align:
			sh = np.array([y_shift_cells(g, i, j, jn, False, dz)
				for i in range(g.nx)])[x]
			b = shift_z(b, align * sh)
		return a, b

	ne0, ne1 = pair(s["ne"])
	nef = 0.5 * (ne0 + ne1)
	L = 0.5 * (np.sqrt(g["g_22"][x, j]) * g["dy"][x, j]
		+ np.sqrt(g["g_22"][x, jn]) * g["dy"][x, jn])[:, None]
	area = (g["J_ylow"][x, jn] * g["dx_ylow"][x, jn] * dz
		/ np.sqrt(g["g_22_ylow"][x, jn]))[:, None]

	te0, te1 = pair(s["te"])
	ti0, ti1 = pair(s["ti"])
	nu_e0, nu_i0 = collision_frequencies(ne0, te0, ti0, ai)
	nu_e1, nu_i1 = collision_frequencies(ne1, te1, ti1, ai)

	out = {}
	for sp, T0, T1, nu0, nu1, A, coef, v in (
			("e", te0, te1, nu_e0, nu_e1, ae, pars["kappa_coefficient_e"],
				s["ve"]),
			("i", ti0, ti1, nu_i0, nu_i1, ai, pars["kappa_coefficient_i"],
				s["vi"])):
		# kappa = c P tau / m at each cell, averaged onto the face.
		m = A * MP
		k0 = coef * np.maximum(ne0 * T0 * QE, 0) / soft_floor(nu0, 1e-10) / m
		k1 = coef * np.maximum(ne1 * T1 * QE, 0) / soft_floor(nu1, 1e-10) / m
		q = -0.5 * (k0 + k1) * (T1 - T0) * QE / L
		v0, v1 = pair(v)
		vf = 0.5 * (v0 + v1)
		Tf = 0.5 * (T0 + T1) * QE
		e = 2.5 * nef * Tf * vf + q
		if sp == "i":
			e = e + 0.5 * m * nef * vf ** 3
		out["q_" + sp] = np.sum(q * area) * info["zperiod"]
		out["qlocal_" + sp] = q
		out["E_" + sp] = np.sum(e * area) * info["zperiod"]
		if sp == "i":
			out["particles"] = np.sum(nef * vf * area) * info["zperiod"]
	return out


def check_balance(g, info, snaps, dz, align):
	"""
	Averages the face fluxes next to both targets over the snapshots and
	compares the totals with the input sources.

	@param g The grid.
	@param info The case's JSON contents.
	@param snaps Snapshot paths.
	@param dz Toroidal cell width (rad).
	@param align Field-alignment sign from check_alignment.
	"""
	print(f"\n== 5. Particle and energy balance ({len(snaps)} snapshots) ==")
	pars = info["pars"]
	norm = parse_norm(info)
	ai, _ = parse_species(pars["ion"])
	ae, _ = parse_species(pars["electron"])
	keys = ("particles", "q_e", "q_i", "E_e", "E_i")
	tot = {a: {k: 0.0 for k in keys} for a in (align, 0)}
	for path in snaps:
		s = load_snapshot(path, norm, ai, ae)
		for a in (align, 0):
			inner = face_fluxes(g, s, info, 0, a, dz)
			outer = face_fluxes(g, s, info, g.ny - 2, a, dz)
			for k in keys:
				# Out of the domain: -y at the inner target, +y at the outer.
				tot[a][k] += (outer[k] - inner[k]) / len(snaps)

	src_p = pars["ne_flux"]
	src_e = 2 * pars["heating_W"]
	print(f"  sources: {src_p:.3g} particles/s, {src_e / 1e3:.0f} kW "
		f"({pars['heating_W'] / 1e3:.0f} kW per species)")
	for a, label in ((align, "along the field line"), (0, "at fixed z")):
		t = tot[a]
		e = t["E_e"] + t["E_i"]
		print(f"  gradients {label}:")
		print(f"    particles to the targets {t['particles']:.3g} /s "
			f"({100 * t['particles'] / src_p:.0f}% of the source)")
		print(f"    energy to the targets {e / 1e3:.1f} kW "
			f"({100 * e / src_e:.0f}% of the heating): electrons "
			f"{t['E_e'] / 1e3:.1f} kW (conducted {t['q_e'] / 1e3:.1f}), "
			f"ions {t['E_i'] / 1e3:.1f} kW (conducted {t['q_i'] / 1e3:.1f})")
	print("  (Some particles and energy also leave through the radial "
		"boundaries, so the totals can be somewhat under 100%.)")

	# Local heat fluxes upstream, where the field line moves furthest in z
	# between planes: compare q at fixed z with q along the field line,
	# cell by cell in the SOL.
	isol = slice(g.ixseps1 - MXG, g.nx - 2 * MXG)
	jm = g.j11 + 1 + np.argmax(g["Rxy"][g.ixseps1, g.j11 + 1:g.j22])
	s = load_snapshot(snaps[-1], norm, ai, ae)
	fa = face_fluxes(g, s, info, jm, align, dz)
	f0 = face_fluxes(g, s, info, jm, 0, dz)
	sh = y_shift_cells(g, g.ixseps1, jm, jm + 1, False, dz)
	print(f"  outboard midplane face (y {jm} -> {jm + 1}, field line moves "
		f"{sh:.0f} z cells), SOL cells, last snapshot:")
	for sp in ("e", "i"):
		qa = fa["qlocal_" + sp][isol]
		q0 = f0["qlocal_" + sp][isol]
		rel = np.sqrt(np.mean((q0 - qa) ** 2) / np.mean(qa ** 2))
		cc = np.corrcoef(qa.ravel(), q0.ravel())[0, 1]
		print(f"    q_{sp}: rms(fixed z - aligned) / rms(aligned) "
			f"{rel:.2f}, correlation {cc:+.2f}")


def main():
	parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
	parser.add_argument("--root", required=True,
		help="Root of the JAX-Hermes3 data (holds inputs/ and data/)")
	parser.add_argument("--case", required=True,
		help="Case name (the JSON file is inputs/<case>.json)")
	parser.add_argument("--grid", required=True,
		help="Grid file of the case, relative to --root or absolute")
	parser.add_argument("--grid-other", default="",
		help="Grid for the opposite field direction, for the grad-B drift "
		"test (empty to skip)")
	parser.add_argument("--nsnap", type=int, default=10,
		help="Number of snapshots, spread over the last --span-ms, to "
		"average the balances over")
	parser.add_argument("--span-ms", type=float, default=2.0)
	args = parser.parse_args()

	g = Grid(os.path.join(args.root, args.grid))
	info, snaps = load_case(args.root, args.case)
	print(f"Case {args.case}: {len(snaps)} snapshots, grid {args.grid} "
		f"({g.nx} x {g.ny}, ixseps1 {g.ixseps1}, jyseps1_1 {g.j11}, "
		f"jyseps2_2 {g.j22}), zperiod {info['zperiod']}")

	check_corners(g)
	check_handedness(g, args.grid)
	if args.grid_other:
		check_handedness(Grid(os.path.join(args.root, args.grid_other)),
			args.grid_other)

	norm = parse_norm(info)
	ai, _ = parse_species(info["pars"]["ion"])
	ae, _ = parse_species(info["pars"]["electron"])
	last = load_snapshot(snaps[-1], norm, ai, ae)
	nz = last["ne"].shape[2]
	dz = 2 * np.pi / (info["zperiod"] * nz)
	align = check_alignment(g, last["ne"], dz)
	check_flow(g, last, ai)

	# Snapshots over the last span_ms, evenly spaced.
	dt_ms = info["time_unit_seconds"] * 1e3
	t_last = float(np.load(snaps[-1])["t"])
	times = np.array([int(os.path.basename(p)[4:-4]) for p in snaps])
	pick = times >= times[-1] - args.span_ms / dt_ms
	chosen = np.array(snaps)[pick]
	idx = np.linspace(0, len(chosen) - 1, args.nsnap).round().astype(int)
	print(f"\n(last snapshot t = {t_last:.0f} t.u. = "
		f"{t_last * dt_ms:.3f} ms)")
	check_balance(g, info, list(chosen[np.unique(idx)]), dz, align)


if __name__ == "__main__":
	main()
