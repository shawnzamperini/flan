"""
Checks the conventions of a HERMES run (JAX-Hermes3 or BOUT++ Hermes-3)
before it is used as a Flan background.

A JAX-Hermes3 run is a hypnotoad grid (netCDF) plus one .npz snapshot per
time with Ne, Pe, Pi, NVi, NVe, phi and Vort on the (x, y, z) cell centers,
all normalized (Nnorm, Tnorm, Bnorm). The case's JSON file gives the
normalization, the species masses, the conduction coefficients, the
sources and zperiod.

A BOUT++ Hermes-3 run is the same kind of grid plus one dump file per MPI
process (BOUT.dmp.N.nc), each holding its block of cells (with guard
cells) at every output time, and the normalization. The species masses
and conduction coefficients come from BOUT.settings, and the sources from
the dumped source diagnostics (Se_src, Pe_src, Pi_src). The blocks are
assembled into global fields as Flan's reader does: interior cells, plus
x guard cells at the domain's x boundaries, no y guard cells.

Flan's HERMES reader relies on several conventions that the files don't
state, so this script tests each one on the data:

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

For BOUT++ runs there is also:

6. Dump checks. (a) Every global cell is filled by exactly one block, and
   field-aligned differences across block boundaries look like those
   inside blocks. (b) In an electromagnetic run Hermes-3 evolves the
   canonical momentum, so the dumped NV = AA N V + Z N Apar; this fits
   NV - AA N V against N Apar. (c) The kappa_par this script calculates
   against the dumped kappa_par_e/i, and its conducted power through a
   few y faces against the dumped efe/efi_cond_ylow. (d) The inductive
   parallel electric field -dApar/dt (from neighboring output times)
   against the electrostatic -grad_par phi, which is all Flan uses.

Usage:

	python check_hermes_conventions.py --root /path/to/data \\
		--case my_case --grid mesh/reversed.nc

	python check_hermes_conventions.py --format bout \\
		--run-dir /path/to/run [--grid grid.nc]

The forward-field grid is used for the grad-B drift test as well if
--grid-other is given (e.g. mesh/forward.nc). For BOUT++ runs the grid
defaults to the [mesh] file in BOUT.settings, and --cache saves the
assembled fields to a .npz file (and reads them back if it exists), since
assembling them from many dump files takes minutes.
"""

import argparse
import ast
import glob
import json
import operator
import os
import re
import time

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


class JaxCase:
	"""
	A JAX-Hermes3 case: the parameters the checks need, from the case's
	JSON file, and its snapshots, read one at a time.
	"""

	def __init__(self, root, case):
		"""
		@param root Root of the JAX-Hermes3 data (holds inputs/ and data/).
		@param case Case name (the JSON file is inputs/<case>.json).
		"""
		self.info, self.paths = load_case(root, case)
		pars = self.info["pars"]
		self.norm = parse_norm(self.info)
		self.ai, _ = parse_species(pars["ion"])
		self.ae, _ = parse_species(pars["electron"])
		self.kappa_e = pars["kappa_coefficient_e"]
		self.kappa_i = pars["kappa_coefficient_i"]
		self.zperiod = self.info["zperiod"]
		self.time_unit_s = self.info["time_unit_seconds"]
		# Sources from the input: heating_W goes to each species.
		self.src_particles = pars["ne_flux"]
		self.src_power = {"e": pars["heating_W"], "i": pars["heating_W"]}

	def snapshot(self, path):
		"""
		@param path Path of a snapshot (from self.paths).
		@return The snapshot in SI units (see load_snapshot).
		"""
		return load_snapshot(path, self.norm, self.ai, self.ae)


# ---------------------------------------------------------------------------
# BOUT++ Hermes-3 dumps
# ---------------------------------------------------------------------------

# Operators allowed in BOUT.settings values evaluated as numbers.
_OPS = {ast.Add: operator.add, ast.Sub: operator.sub, ast.Mult: operator.mul,
	ast.Div: operator.truediv, ast.Pow: operator.pow}


def eval_number(text):
	"""
	Evaluates a BOUT++ option value that is plain arithmetic on numbers,
	such as "1 / 1836" (+ - * / and ^ only, no functions or names).

	@param text The option's value.
	@return The value as a float.
	@throws ValueError If it isn't plain arithmetic on numbers.
	"""
	def ev(n):
		if isinstance(n, ast.Constant) and isinstance(n.value, (int, float)):
			return float(n.value)
		if isinstance(n, ast.UnaryOp) and isinstance(n.op, (ast.UAdd,
				ast.USub)):
			v = ev(n.operand)
			return v if isinstance(n.op, ast.UAdd) else -v
		if isinstance(n, ast.BinOp) and type(n.op) in _OPS:
			return _OPS[type(n.op)](ev(n.left), ev(n.right))
		raise ValueError(f"not a number: {text}")
	return ev(ast.parse(text.replace("^", "**"), mode="eval").body)


def parse_settings(path):
	"""
	Reads a BOUT.settings (or BOUT.inp) file: "[section]" headers and
	"key = value" lines, with "#" comments.

	@param path Path to the file.
	@return Dict of sections ("" for the top level), each a dict of
	        key -> value string.
	"""
	out = {"": {}}
	sec = ""
	with open(path) as f:
		for line in f:
			line = line.split("#", 1)[0].strip()
			if not line:
				continue
			if line.startswith("[") and line.endswith("]"):
				sec = line[1:-1].strip()
				out.setdefault(sec, {})
			elif "=" in line:
				k, v = line.split("=", 1)
				out[sec][k.strip()] = v.strip()
	return out


class BoutRun:
	"""
	A BOUT++ Hermes-3 run: the process layout and normalization from the
	first dump file, the species parameters from BOUT.settings, and the
	assembly of global (x, y, z) fields from all the dump files.
	"""

	def __init__(self, run_dir, ion="i"):
		"""
		@param run_dir Directory with BOUT.dmp.*.nc and BOUT.settings.
		@param ion Name of the ion species in the Hermes-3 input.
		"""
		self.dir = run_dir
		self.ion = ion
		with netCDF4.Dataset(os.path.join(run_dir, "BOUT.dmp.0.nc")) as d:
			def sc(k):
				return d[k][...].item()
			(self.nxpe, self.nype, self.nzpe, self.mxsub, self.mysub,
				self.mxg, self.myg, self.nz) = (int(sc(k)) for k in
				("NXPE", "NYPE", "NZPE", "MXSUB", "MYSUB", "MXG", "MYG",
				"MZ"))
			self.zperiod = int(round(sc("zperiod")))
			self.norm = {k: float(sc(k)) for k in ("Nnorm", "Tnorm",
				"Bnorm", "Cs0", "Omega_ci", "rho_s0")}
			self.seps = {k: int(sc(k)) for k in ("nx", "ny", "ixseps1",
				"jyseps1_1", "jyseps2_1", "jyseps1_2", "jyseps2_2")}
			self.t = np.array(d["t_array"][:], dtype=float)
			self.vars = set(d.variables)
			# SI value = dumped value * conversion
			self.conv = {k: float(d[k].getncattr("conversion"))
				for k in d.variables if "conversion" in d[k].ncattrs()}
		if self.nzpe != 1:
			raise SystemExit(f"NZPE = {self.nzpe}: z split over processes "
				"isn't handled.")
		self.nx = self.nxpe * self.mxsub + 2 * self.mxg
		self.ny = self.nype * self.mysub
		self.paths = [os.path.join(run_dir, f"BOUT.dmp.{p}.nc")
			for p in range(self.nxpe * self.nype)]
		missing = [p for p in self.paths if not os.path.exists(p)]
		if missing:
			raise SystemExit(f"{len(missing)} dump files are missing, e.g. "
				f"{missing[0]}.")
		self.time_unit_s = 1 / self.norm["Omega_ci"]

		self.settings = parse_settings(os.path.join(run_dir, "BOUT.settings"))
		if ion not in self.settings or "e" not in self.settings:
			raise SystemExit(f"BOUT.settings has no [{ion}] or [e] section; "
				"set --ion-species to the ion species' name.")
		si, se = self.settings[ion], self.settings["e"]
		self.ai = eval_number(si["AA"])
		self.zi = eval_number(si.get("charge", "1"))
		self.ae = eval_number(se["AA"])
		# Hermes-3's defaults: 3.16 / sqrt(2) for electrons, 3.9 for ions.
		self.kappa_e = eval_number(se.get("kappa_coefficient",
			str(3.16 / np.sqrt(2))))
		self.kappa_i = eval_number(si.get("kappa_coefficient", "3.9"))
		self.limit_e = eval_number(se.get("kappa_limit_alpha", "-1"))
		self.limit_i = eval_number(si.get("kappa_limit_alpha", "-1"))
		self.collisions = {sp: self.settings[sp].get(
			"conduction_collisions_mode", "multispecies") for sp in ("e", ion)}
		comps = self.settings.get("hermes", {}).get("components", "")
		self.em = "electromagnetic" in comps or "Apar" in self.vars
		self.grid_name = self.settings.get("mesh", {}).get("file", "") \
			.strip("\"'")
		self.data = {}

	def check_grid(self, g):
		"""
		Stops if the dumps' sizes or separatrix indices differ from the
		grid's.

		@param g The grid.
		"""
		want = {"nx": g.nx, "ny": g.ny, "ixseps1": g.ixseps1,
			"jyseps1_1": g.j11, "jyseps2_2": g.j22}
		have = dict(self.seps, nx=self.nx, ny=self.ny)
		bad = {k: (have[k], v) for k, v in want.items() if have[k] != v}
		if self.seps["nx"] != self.nx or self.seps["ny"] != self.ny:
			bad["layout"] = ((self.nx, self.ny),
				(self.seps["nx"], self.seps["ny"]))
		if bad:
			raise SystemExit("The dumps and the grid disagree (dumps, grid): "
				+ ", ".join(f"{k} {a} vs {b}" for k, (a, b) in bad.items()))

	def read(self, requests):
		"""
		Assembles the requested variables at the requested output times
		from all the dump files, opening each file once.

		From each block the interior cells are taken, plus the x guard
		cells where the block is at the domain's x boundary; y guard cells
		are dropped (the grid has none). Stops if a global cell isn't
		filled exactly once, or a file's times differ from the first's.

		@param requests Dict of variable name -> list of time indices.
		@return Nothing; self.data[name] = (time indices, array (nt, nx,
		        ny, nz)) of the dumped (normalized) values.
		"""
		reqs = {v: sorted(set(fr)) for v, fr in requests.items()}
		for v in reqs:
			if v not in self.vars:
				raise SystemExit(f"The dumps have no variable {v}.")
		out = {v: np.full((len(fr), self.nx, self.ny, self.nz), np.nan)
			for v, fr in reqs.items()}
		count = np.zeros((self.nx, self.ny), dtype=int)
		t0 = time.time()
		for n, path in enumerate(self.paths):
			with netCDF4.Dataset(path) as d:
				d.set_auto_mask(False)
				if not np.array_equal(d["t_array"][:], self.t):
					raise SystemExit(f"{path} has different output times "
						"from BOUT.dmp.0.nc.")
				px, py = int(d["PE_XIND"][...]), int(d["PE_YIND"][...])
				x0 = 0 if px == 0 else self.mxg
				x1 = self.mxsub + self.mxg + (self.mxg if px == self.nxpe - 1
					else 0)
				gx = px * self.mxsub + x0
				gy = py * self.mysub
				count[gx:gx + x1 - x0, gy:gy + self.mysub] += 1
				for v, fr in reqs.items():
					out[v][:, gx:gx + x1 - x0, gy:gy + self.mysub] = \
						d[v][fr, x0:x1, self.myg:self.myg + self.mysub, :]
			if (n + 1) % 128 == 0:
				print(f"  read {n + 1} of {len(self.paths)} dump files "
					f"({time.time() - t0:.0f} s)", flush=True)
		if (count != 1).any():
			raise SystemExit(f"{(count == 0).sum()} global cells weren't "
				f"filled and {(count > 1).sum()} were filled more than once.")
		self.data = {v: (reqs[v], out[v]) for v in reqs}

	def save(self, path):
		"""
		Saves the assembled fields to a .npz file (see load).

		@param path Path of the file.
		"""
		np.savez(path, **{f"{v}": a for v, (_, a) in self.data.items()},
			**{f"frames:{v}": np.array(fr) for v, (fr, _) in
			self.data.items()})

	def load(self, path, requests):
		"""
		Reads fields saved by save, if they hold everything requested.

		@param path Path of the file.
		@param requests As for read.
		@return True if the file held everything requested.
		"""
		if not os.path.exists(path):
			return False
		z = np.load(path)
		data = {}
		for v, fr in requests.items():
			if f"frames:{v}" not in z or not set(fr) <= set(z[f"frames:{v}"]):
				return False
			data[v] = (list(z[f"frames:{v}"]), z[v])
		self.data = data
		return True

	def get(self, name, frame):
		"""
		@param name Variable name.
		@param frame Output time index.
		@return The assembled, dumped (normalized) field (nx, ny, nz).
		"""
		frames, a = self.data[name]
		return a[frames.index(frame)]

	def has(self, name):
		"""@return Whether the dumps have the variable."""
		return name in self.vars

	def flow_vars(self, sp):
		"""
		Variables needed for a species' parallel flow: V if it was dumped,
		else NV, and Apar in an electromagnetic run.

		@param sp Species name.
		@return List of variable names.
		"""
		if self.has("V" + sp):
			return ["V" + sp]
		return ["NV" + sp] + (["Apar"] if self.em else [])

	def flow(self, sp, frame, aa, z, nf):
		"""
		A species' parallel flow, normalized (units of Cs0): the dumped V,
		or (NV - Z N Apar) / (AA N), NV being the canonical momentum in an
		electromagnetic run.

		@param sp Species name.
		@param frame Output time index.
		@param aa, z The species' mass number and charge.
		@param nf Floored density (normalized).
		@return The flow (nx, ny, nz).
		"""
		if self.has("V" + sp):
			return self.get("V" + sp, frame)
		nv = self.get("NV" + sp, frame)
		if self.em:
			nv = nv - z * self.get("Ne", frame) * self.get("Apar", frame)
		return nv / (aa * nf)

	def snapshot_vars(self):
		"""@return Variables needed by snapshot."""
		return ["Ne", "Pe", "P" + self.ion] + self.flow_vars(self.ion) \
			+ self.flow_vars("e")

	def snapshot(self, frame):
		"""
		One output time in SI units, like load_snapshot for JAX-Hermes3.
		Needs snapshot_vars() read at that time.

		@param frame Output time index.
		@return Dict with t, ne (m^-3), te, ti (eV), vi, ve (m/s).
		"""
		n = self.get("Ne", frame)
		nf = soft_floor(n, 1e-5)
		cs0 = self.norm["Cs0"]
		return {
			"t": float(self.t[frame]),
			"ne": n * self.norm["Nnorm"],
			"te": self.get("Pe", frame) / nf * self.norm["Tnorm"],
			"ti": self.get("P" + self.ion, frame) / nf * self.norm["Tnorm"],
			"vi": self.flow(self.ion, frame, self.ai, self.zi, nf) * cs0,
			"ve": self.flow("e", frame, self.ae, -1.0, nf) * cs0,
		}

	def set_sources(self, g, frames, dz):
		"""
		Total particle and heating sources (averaged over the given
		times) from the dumped external sources: the integral of Se_src,
		and of 3/2 P_src for each species, over the physical cells and the
		whole torus, with the cell volume J dx dy dz from the grid.

		@param g The grid.
		@param frames Output time indices (Se_src, Pe_src and Pi_src read).
		@param dz Toroidal cell width (rad).
		"""
		x = slice(MXG, g.nx - MXG)
		vol = (g["J"] * g["dx"] * g["dy"])[x, :, None] * dz * self.zperiod
		def total(name):
			return np.mean([np.sum(self.get(name, f)[x] * vol)
				for f in frames]) * self.conv[name]
		self.src_particles = total("Se_src")
		self.src_power = {"e": 1.5 * total("Pe_src"),
			"i": 1.5 * total(f"P{self.ion}_src")}


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


def shift_spectral(f, s):
	"""
	Samples a periodic array at k + s along its last axis, through its
	trigonometric interpolant (as BOUT++'s shifted metric and Flan's
	aligned derivatives do), so resolved structure isn't smoothed.

	@param f Array (..., nz).
	@param s Shift in cells (scalar, or broadcastable to f[..., 0]).
	@return f at k + s for each k.
	"""
	nz = f.shape[-1]
	m = np.fft.fftfreq(nz) * nz
	ph = np.exp(2j * np.pi * m * np.asarray(s)[..., None] / nz)
	return np.fft.ifft(np.fft.fft(f, axis=-1) * ph, axis=-1).real


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


def face_fluxes(g, s, case, j, align, dz):
	"""
	Particle flux and parallel energy flux (W) through the y face between
	cells j and j + 1, summed over x (physical cells) and the whole torus,
	positive toward +y.

	The face area normal to b is J dx dz / sqrt(g_22) at the face (the
	_ylow fields of cell j + 1), times zperiod for the full torus. Face
	values are the average of the two cells; the temperature gradient is
	(T(j + 1) - T(j)) / L, with L = (sqrt(g_22) dy) averaged over the two
	cells. With align = +/-1, the j + 1 values are taken along the field
	line (z +/- the zShift change), shifted spectrally as BOUT++ and Flan
	do (a linear shift smooths the plane and underestimates gradients);
	with align = 0 at fixed z.

	The energy flux per species is 5/2 n T V + q, plus 1/2 m_i n V^3 for
	the ions, with q = -kappa grad_par T and kappa from
	braginskii_conduction.

	@param g The grid.
	@param s A snapshot in SI units.
	@param case The case (JaxCase or BoutRun): masses, conduction
	       coefficients and zperiod.
	@param j y index of the cell before the face.
	@param align +1, -1 or 0, as above.
	@param dz Toroidal cell width (rad).
	@return Dict with particles (s^-1), and conducted and total energy
	        fluxes (W) for each species, plus the local conducted heat
	        fluxes (qlocal_e, qlocal_i, W/m^2) and face areas (area, m^2)
	        of the physical x cells.
	"""
	ai, ae = case.ai, case.ae
	x = slice(MXG, g.nx - MXG)
	jn = j + 1

	def pair(f):
		a = f[x, j]
		b = f[x, jn]
		if align:
			sh = np.array([y_shift_cells(g, i, j, jn, False, dz)
				for i in range(g.nx)])[x]
			b = shift_spectral(b, align * sh)
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
			("e", te0, te1, nu_e0, nu_e1, ae, case.kappa_e, s["ve"]),
			("i", ti0, ti1, nu_i0, nu_i1, ai, case.kappa_i, s["vi"])):
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
		out["q_" + sp] = np.sum(q * area) * case.zperiod
		out["qlocal_" + sp] = q
		out["E_" + sp] = np.sum(e * area) * case.zperiod
		if sp == "i":
			out["particles"] = np.sum(nef * vf * area) * case.zperiod
	out["area"] = area
	return out


def check_balance(g, case, snaps, dz, align):
	"""
	Averages the face fluxes next to both targets over the snapshots and
	compares the totals with the sources.

	@param g The grid.
	@param case The case (JaxCase or BoutRun), with snapshot(key) and the
	       sources src_particles (s^-1) and src_power (W per species).
	@param snaps Snapshot keys for case.snapshot, in time order.
	@param dz Toroidal cell width (rad).
	@param align Field-alignment sign from check_alignment.
	"""
	print(f"\n== 5. Particle and energy balance ({len(snaps)} snapshots) ==")
	keys = ("particles", "q_e", "q_i", "E_e", "E_i")
	tot = {a: {k: 0.0 for k in keys} for a in (align, 0)}
	for key in snaps:
		s = case.snapshot(key)
		for a in (align, 0):
			inner = face_fluxes(g, s, case, 0, a, dz)
			outer = face_fluxes(g, s, case, g.ny - 2, a, dz)
			for k in keys:
				# Out of the domain: -y at the inner target, +y at the outer.
				tot[a][k] += (outer[k] - inner[k]) / len(snaps)

	src_p = case.src_particles
	src_e = case.src_power["e"] + case.src_power["i"]
	print(f"  sources: {src_p:.3g} particles/s, {src_e / 1e3:.0f} kW "
		f"(electrons {case.src_power['e'] / 1e3:.0f} kW, ions "
		f"{case.src_power['i'] / 1e3:.0f} kW)")
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
		"boundaries, so the totals can be somewhat under 100%; they can "
		"be over 100% if the plasma isn't steady, e.g. its density is "
		"falling.)")

	# Local heat fluxes upstream, where the field line moves furthest in z
	# between planes: compare q at fixed z with q along the field line,
	# cell by cell in the SOL.
	isol = slice(g.ixseps1 - MXG, g.nx - 2 * MXG)
	jm = g.j11 + 1 + np.argmax(g["Rxy"][g.ixseps1, g.j11 + 1:g.j22])
	s = case.snapshot(snaps[-1])
	fa = face_fluxes(g, s, case, jm, align, dz)
	f0 = face_fluxes(g, s, case, jm, 0, dz)
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


# ---------------------------------------------------------------------------
# 6. BOUT++ dump checks
# ---------------------------------------------------------------------------

def aligned_neighbor(g, f, j, jn, cut, dz, align):
	"""
	Plane jn of f, shifted onto the field lines through the cells of
	plane j.

	@param g The grid.
	@param f Field (nx, ny, nz).
	@param j, jn y indices (jn = j + 1 or j - 1, or across a cut).
	@param cut Whether going from j to jn crosses the core branch cut.
	@param dz Toroidal cell width (rad).
	@param align Field-alignment sign from check_alignment.
	@return Array (nx, nz).
	"""
	if jn > j or cut:
		sh = np.array([y_shift_cells(g, i, j, jn, cut, dz)
			for i in range(g.nx)])
	else:
		sh = -np.array([y_shift_cells(g, i, jn, j, cut, dz)
			for i in range(g.nx)])
	return shift_spectral(f[:, jn], align * sh)


def check_assembly(g, run, ne, dz, align):
	"""
	Compares field-aligned differences of ln(Ne) across the process
	block boundaries with those inside the blocks, in y and in x. With
	the blocks put in the right places the two are alike; a misplaced
	block shows up as much larger jumps at its boundaries.

	@param g The grid.
	@param run The BoutRun (filled exactly once is checked when reading).
	@param ne Electron density (nx, ny, nz), any units.
	@param dz Toroidal cell width (rad).
	@param align Field-alignment sign from check_alignment.
	"""
	print("\n== 6a. Assembly of the dump blocks ==")
	print(f"  {run.nxpe} x {run.nype} blocks of {run.mxsub} x {run.mysub} "
		f"cells (+{run.mxg} x, {run.myg} y guards): every one of the "
		f"{g.nx} x {g.ny} cells filled once")
	lnn = np.log(np.maximum(ne, 1e-30))
	x = slice(MXG, g.nx - MXG)
	jumps = {True: [], False: []}
	for j in range(g.ny - 1):
		# Faces with the same neighbor across them at every x (no cut).
		if j in (g.j11, g.j22):
			continue
		d = np.abs(aligned_neighbor(g, lnn, j, j + 1, False, dz, align)
			- lnn[:, j])[x].mean()
		jumps[(j + 1) % run.mysub == 0].append(d)
	print(f"  y faces, mean |d ln Ne| along the field line: block "
		f"boundaries {np.mean(jumps[True]):.4f}, inside blocks "
		f"{np.mean(jumps[False]):.4f}")
	dx = np.abs(np.diff(lnn[MXG:g.nx - MXG], axis=0)).mean(axis=(1, 2))
	at = np.array([(i + 1) % run.mxsub == 0 for i in range(len(dx))])
	print(f"  x faces, mean |d ln Ne|: block boundaries {dx[at].mean():.4f}, "
		f"inside blocks {dx[~at].mean():.4f}")


def check_canonical(run, frame):
	"""
	Fits NV - AA N V against N Apar for each species, to find whether the
	dumped NV is the canonical momentum (slope Z) or the plain momentum
	(no Apar, or slope 0), all in the dumps' normalized units.

	@param run The BoutRun, with Ne, Apar and each species' NV and V read.
	@param frame Output time index.
	"""
	print("\n== 6b. Momentum and flow ==")
	x = slice(run.mxg, run.nx - run.mxg)
	n = run.get("Ne", frame)[x]
	for sp, aa, z in (("e", run.ae, -1.0), (run.ion, run.ai, run.zi)):
		if not (run.has("NV" + sp) and run.has("V" + sp)):
			print(f"  {sp}: NV{sp} or V{sp} not dumped, can't check")
			continue
		r = run.get("NV" + sp, frame)[x] - aa * n / abs(z) \
			* run.get("V" + sp, frame)[x]
		rel = np.sqrt(np.mean(r ** 2)) / np.sqrt(np.mean(
			run.get("NV" + sp, frame)[x] ** 2))
		if not run.has("Apar"):
			print(f"  {sp}: no Apar; rms(NV - AA N V) / rms(NV) = {rel:.2e}")
			continue
		a = (n * run.get("Apar", frame)[x]).ravel()
		slope, icpt = np.polyfit(a, r.ravel(), 1)
		res = r.ravel() - (slope * a + icpt)
		print(f"  {sp}: NV - AA N V = {slope:+.4f} N Apar {icpt:+.1e}, "
			f"residual {np.std(res) / max(np.std(r), 1e-300):.1e} of its "
			f"rms; NV differs from AA N V by {rel:.1%} rms (Z = {z:+g})")
	print("  -> the flow has to come from V (or NV - Z N Apar), not "
		"NV / (AA N)" if run.has("Apar") else "")


def check_conduction(g, run, frames, dz, align):
	"""
	Compares this script's conduction with the run's: kappa_par per cell
	against the dumped kappa_par_e/i, and the conducted power through a
	few y faces against the dumped efe/efi_cond_ylow (the power through
	each cell's lower y face, positive toward +y).

	The face powers are compared summed over z for each x, which doesn't
	depend on which toroidal frame the face values are in, and over all
	x. Hermes-3 interpolates to the face in its own way, so a few percent
	of difference is expected; a sign or a large factor isn't.

	@param g The grid.
	@param run The BoutRun, with the snapshot variables, kappa_par_e/i and
	       efe/efi_cond_ylow read at frames.
	@param frames Output time indices to average over.
	@param dz Toroidal cell width (rad).
	@param align Field-alignment sign from check_alignment.
	"""
	print(f"\n== 6c. Conduction against the run's diagnostics "
		f"({len(frames)} times) ==")
	print(f"  BOUT.settings: kappa_coefficient e {run.kappa_e:.5g}, "
		f"{run.ion} {run.kappa_i:.5g}; kappa_limit_alpha e {run.limit_e:g}, "
		f"{run.ion} {run.limit_i:g}; collisions {run.collisions}")
	if run.limit_e > 0 or run.limit_i > 0:
		print("  !! A flux limiter is on; this script's kappa has none, so "
			"expect differences where it acts.")
	x = slice(MXG, g.nx - MXG)
	ion = run.ion
	s = run.snapshot(frames[-1])
	nu_e, nu_i = collision_frequencies(s["ne"], s["te"], s["ti"], run.ai)
	for sp, nu, T, coef, aa in (("e", nu_e, s["te"], run.kappa_e, run.ae),
			(ion, nu_i, s["ti"], run.kappa_i, run.ai)):
		name = f"kappa_par_{sp}"
		if not run.has(name):
			print(f"  {name} not dumped")
			continue
		# W/m/eV: c n T tau / m, times e for T in eV.
		mine = coef * s["ne"] * T * QE / soft_floor(nu, 1e-10) / (aa * MP) \
			* QE
		theirs = run.get(name, frames[-1]) * run.conv[name]
		r = (mine / theirs)[x]
		print(f"  {name}: this script / dumped, median {np.median(r):.4f}, "
			f"1-99% [{np.percentile(r, 1):.4f}, {np.percentile(r, 99):.4f}]")

	jm = g.j11 + 1 + np.argmax(g["Rxy"][g.ixseps1, g.j11 + 1:g.j22])
	for j, label in ((0, "next to the inner target"),
			(jm, "outboard midplane"), (g.ny - 2, "next to the outer target")):
		tot = {"e": [0.0, 0.0], "i": [0.0, 0.0]}
		rows = {"e": [], "i": []}
		for f in frames:
			s = run.snapshot(f)
			ff = face_fluxes(g, s, run, j, align, dz)
			for sp, dname in (("e", "efe_cond_ylow"), ("i", "efi_cond_ylow")):
				mine = (ff["qlocal_" + sp] * ff["area"]).sum(axis=1)
				theirs = (run.get(dname, f)[:, j + 1] * run.conv[dname]) \
					[x].sum(axis=1)
				tot[sp][0] += mine.sum() / len(frames)
				tot[sp][1] += theirs.sum() / len(frames)
				rows[sp].append((mine, theirs))
		print(f"  y face {j} -> {j + 1} ({label}), conducted power over the "
			"wedge:")
		for sp in ("e", "i"):
			mine = np.concatenate([m for m, _ in rows[sp]])
			theirs = np.concatenate([t for _, t in rows[sp]])
			cc = np.corrcoef(mine, theirs)[0, 1]
			print(f"    {sp}: this script {tot[sp][0] / 1e3:+9.2f} kW, dumped "
				f"{tot[sp][1] / 1e3:+9.2f} kW; per-x correlation {cc:+.3f}")


def check_inductive(g, run, centers, dz, align):
	"""
	Compares the inductive parallel electric field -dApar/dt with the
	electrostatic one, -grad_par phi, in the core and the SOL.

	dApar/dt is a centered difference over neighboring output times; the
	same over twice the interval shows whether the output times resolve
	it (if the two differ a lot, Apar changes faster than the output
	interval and dApar/dt here is only a rough size). grad_par phi is a
	centered difference along the field line,
	(phi(j + 1) - phi(j - 1)) / (distance between the cell centers),
	with the neighboring planes shifted onto the field line. The target
	cells, and core cells next to the branch cut, are left out (in the
	SOL, y runs from target to target without a cut). Both the full fields and
	their toroidal averages are compared: the averages are what push
	impurities the same way over many turbulent structures.

	@param g The grid.
	@param run The BoutRun, with phi read at centers and Apar at centers
	       +/- 1 and +/- 2.
	@param centers Output time indices.
	@param dz Toroidal cell width (rad).
	@param align Field-alignment sign from check_alignment.
	"""
	print(f"\n== 6d. Inductive vs electrostatic E_par ({len(centers)} "
		"times) ==")
	if not run.has("Apar"):
		print("  No Apar in the dumps: electrostatic run, nothing to check.")
		return
	L = np.sqrt(g["g_22"]) * g["dy"]
	ca, cp = run.conv["Apar"], run.conv["phi"]
	core_y = set(range(g.j11 + 2, g.j22))
	res = {r: {k: [] for k in ("ind", "es", "ind_avg", "es_avg", "c12")}
		for r in ("core", "SOL")}
	for c in centers:
		dt1 = (run.t[c + 1] - run.t[c - 1]) * run.time_unit_s
		dt2 = (run.t[c + 2] - run.t[c - 2]) * run.time_unit_s
		e1 = -(run.get("Apar", c + 1) - run.get("Apar", c - 1)) * ca / dt1
		e2 = -(run.get("Apar", c + 2) - run.get("Apar", c - 2)) * ca / dt2
		phi = run.get("phi", c) * cp
		for region, xs, ys in (
				("core", range(MXG, g.ixseps1), sorted(core_y)),
				("SOL", range(g.ixseps1, g.nx - MXG), range(1, g.ny - 1))):
			xs = list(xs)
			for j in ys:
				up = aligned_neighbor(g, phi, j, j + 1, False, dz, align)
				dn = aligned_neighbor(g, phi, j, j - 1, False, dz, align)
				dist = 0.5 * L[:, j - 1] + L[:, j] + 0.5 * L[:, j + 1]
				es = -((up - dn) / dist[:, None])[xs]
				ind = e1[xs, j]
				res[region]["es"].append(es)
				res[region]["ind"].append(ind)
				res[region]["es_avg"].append(es.mean(axis=-1))
				res[region]["ind_avg"].append(ind.mean(axis=-1))
				res[region]["c12"].append((ind, e2[xs, j]))

	def rms(lst):
		return np.sqrt(np.mean(np.concatenate([a.ravel() for a in lst]) ** 2))

	for region, r in res.items():
		a = np.concatenate([p.ravel() for p, _ in r["c12"]])
		b = np.concatenate([q.ravel() for _, q in r["c12"]])
		cc = np.corrcoef(a, b)[0, 1]
		print(f"  {region}: rms -dApar/dt {rms(r['ind']):.3g} V/m, "
			f"rms -grad_par phi {rms(r['es']):.3g} V/m (ratio "
			f"{rms(r['ind']) / rms(r['es']):.3f})")
		print(f"  {'':{len(region)}}  toroidal averages: {rms(r['ind_avg']):.3g} "
			f"vs {rms(r['es_avg']):.3g} V/m (ratio "
			f"{rms(r['ind_avg']) / rms(r['es_avg']):.3f}); dApar/dt over 1 "
			f"and 2 output intervals correlate {cc:+.3f}")


def pick_balance(times, time_unit_s, span_ms, nsnap):
	"""
	Indices of nsnap times evenly spread over the last span_ms (fewer if
	there aren't that many), always including the last.

	@param times Output times, in the run's time unit, increasing.
	@param time_unit_s The time unit in seconds.
	@param span_ms Span to spread them over (ms).
	@param nsnap Number of times.
	@return Sorted list of indices into times.
	"""
	times = np.asarray(times)
	chosen = np.flatnonzero(times >= times[-1] - span_ms * 1e-3 / time_unit_s)
	idx = np.linspace(0, len(chosen) - 1, nsnap).round().astype(int)
	return [int(i) for i in chosen[np.unique(idx)]]


def main_jax(args):
	"""
	Runs the checks on a JAX-Hermes3 case.

	@param args Parsed command line.
	"""
	g = Grid(os.path.join(args.root, args.grid))
	case = JaxCase(args.root, args.case)
	snaps = case.paths
	print(f"Case {args.case}: {len(snaps)} snapshots, grid {args.grid} "
		f"({g.nx} x {g.ny}, ixseps1 {g.ixseps1}, jyseps1_1 {g.j11}, "
		f"jyseps2_2 {g.j22}), zperiod {case.zperiod}")

	check_corners(g)
	check_handedness(g, args.grid)
	if args.grid_other:
		check_handedness(Grid(os.path.join(args.root, args.grid_other)),
			args.grid_other)

	last = case.snapshot(snaps[-1])
	nz = last["ne"].shape[2]
	dz = 2 * np.pi / (case.zperiod * nz)
	align = check_alignment(g, last["ne"], dz)
	check_flow(g, last, case.ai)

	# Snapshots over the last span_ms, evenly spaced. The file names are
	# the times in the case's time unit.
	times = np.array([int(os.path.basename(p)[4:-4]) for p in snaps])
	pick = pick_balance(times, case.time_unit_s, args.span_ms, args.nsnap)
	print(f"\n(last snapshot t = {last['t']:.0f} t.u. = "
		f"{last['t'] * case.time_unit_s * 1e3:.3f} ms)")
	check_balance(g, case, [snaps[i] for i in pick], dz, align)


def main_bout(args):
	"""
	Runs the checks on a BOUT++ Hermes-3 run.

	@param args Parsed command line.
	"""
	run = BoutRun(args.run_dir, args.ion_species)
	grid = args.grid or run.grid_name
	if not grid:
		raise SystemExit("No --grid given and no [mesh] file in "
			"BOUT.settings.")
	if not os.path.isabs(grid):
		grid = os.path.join(args.run_dir, grid)
	g = Grid(grid)
	run.check_grid(g)
	nt = len(run.t)
	print(f"Run {args.run_dir}: {run.nxpe * run.nype} dump files, {nt} "
		f"output times (t = {run.t[0]:.0f} to {run.t[-1]:.0f} / Omega_ci = "
		f"{run.t[0] * run.time_unit_s * 1e3:.4f} to "
		f"{run.t[-1] * run.time_unit_s * 1e3:.4f} ms)")
	print(f"  grid {grid} ({g.nx} x {g.ny}, ixseps1 {g.ixseps1}, jyseps1_1 "
		f"{g.j11}, jyseps2_2 {g.j22}), nz {run.nz}, zperiod {run.zperiod}")
	print(f"  Nnorm {run.norm['Nnorm']:.4g} m^-3, Tnorm {run.norm['Tnorm']:g} "
		f"eV, Bnorm {run.norm['Bnorm']:g} T; ion species '{run.ion}' AA "
		f"{run.ai:g}, charge {run.zi:g}; electron AA {run.ae:.4g}; "
		f"electromagnetic: {'yes' if run.em else 'no'}")
	if run.zi != 1:
		print("  !! Ion charge isn't 1; this script and Flan assume ni = ne.")
	print("  Flan options from BOUT.settings:\n"
		f"    hermes_ion_mass_amu = {run.ai:g}\n"
		f"    hermes_kappa_coef_e = {run.kappa_e:.6g}\n"
		f"    hermes_kappa_coef_i = {run.kappa_i:.6g}\n"
		f"    hermes_kappa_limit_alpha = "
		f"{max(run.limit_e, run.limit_i):g}"
		+ ("  (!! differs between species)" if run.limit_e != run.limit_i
			else ""))
	for sp, mode in run.collisions.items():
		if mode != "multispecies":
			print(f"  !! {sp} conduction_collisions_mode = {mode}; Flan and "
				"this script use multispecies.")

	check_corners(g)
	check_handedness(g, os.path.basename(grid))
	if args.grid_other:
		check_handedness(Grid(args.grid_other), args.grid_other)

	dz = 2 * np.pi / (run.zperiod * run.nz)
	bal = pick_balance(run.t, run.time_unit_s, args.span_ms, args.nsnap)
	# Times for dApar/dt: spread over the run, with two output times on
	# each side.
	centers = sorted(set(np.linspace(2, nt - 3, args.ninductive).round()
		.astype(int).tolist())) if nt >= 5 else []
	req = {}
	def want(names, frames):
		for v in names:
			if run.has(v):
				req.setdefault(v, set()).update(frames)
	want(run.snapshot_vars() + ["Se_src", "Pe_src", f"P{run.ion}_src",
		f"NV{run.ion}", "NVe", f"V{run.ion}", "Ve", "Apar"], bal)
	want(["kappa_par_e", f"kappa_par_{run.ion}", "efe_cond_ylow",
		"efi_cond_ylow"], bal)
	want(["phi"], centers)
	want(["Apar"], [c + d for c in centers for d in (-2, -1, 1, 2)])
	req = {v: sorted(fr) for v, fr in req.items()}
	nbytes = sum(len(f) for f in req.values()) * g.nx * g.ny * run.nz * 8
	if args.cache and run.load(args.cache, req):
		print(f"  assembled fields read from {args.cache}")
	else:
		print(f"  assembling {len(req)} variables ({nbytes / 1e9:.1f} GB) "
			f"from {len(run.paths)} dump files...", flush=True)
		run.read(req)
		if args.cache:
			run.save(args.cache)
			print(f"  saved them to {args.cache}")
	run.set_sources(g, bal, dz)

	last = run.snapshot(bal[-1])
	align = check_alignment(g, last["ne"], dz)
	check_flow(g, last, run.ai)
	print(f"\n(last output t = {last['t']:.0f} / Omega_ci = "
		f"{last['t'] * run.time_unit_s * 1e3:.4f} ms)")
	check_balance(g, run, bal, dz, align)

	check_assembly(g, run, last["ne"], dz, align)
	check_canonical(run, bal[-1])
	if all(run.has(v) for v in ("kappa_par_e", "efe_cond_ylow",
			"efi_cond_ylow")):
		check_conduction(g, run, bal, dz, align)
	else:
		print("\n== 6c. No kappa_par or efe/efi_cond_ylow in the dumps "
			"(diagnose = false), so the conduction can't be compared ==")
	if centers:
		check_inductive(g, run, centers, dz, align)


def main():
	parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
	parser.add_argument("--format", choices=("jax", "bout"), default="jax",
		help="jax: JAX-Hermes3 .npz snapshots; bout: BOUT++ Hermes-3 dumps")
	parser.add_argument("--root", default="",
		help="jax: root of the data (holds inputs/ and data/)")
	parser.add_argument("--case", default="",
		help="jax: case name (the JSON file is inputs/<case>.json)")
	parser.add_argument("--run-dir", default="",
		help="bout: directory with BOUT.dmp.*.nc and BOUT.settings")
	parser.add_argument("--ion-species", default="i",
		help="bout: name of the ion species in the Hermes-3 input")
	parser.add_argument("--grid", default="",
		help="Grid file (jax: relative to --root or absolute; bout: "
		"relative to --run-dir or absolute, default the [mesh] file)")
	parser.add_argument("--grid-other", default="",
		help="Grid for the opposite field direction, for the grad-B drift "
		"test (empty to skip)")
	parser.add_argument("--nsnap", type=int, default=10,
		help="Number of snapshots, spread over the last --span-ms, to "
		"average the balances over")
	parser.add_argument("--span-ms", type=float, default=2.0)
	parser.add_argument("--ninductive", type=int, default=3,
		help="bout: number of times to compare dApar/dt with grad_par phi")
	parser.add_argument("--cache", default="",
		help="bout: .npz file to save the assembled fields to, or read "
		"them from if it has them")
	args = parser.parse_args()

	if args.format == "jax":
		if not (args.root and args.case and args.grid):
			parser.error("--format jax needs --root, --case and --grid")
		main_jax(args)
	else:
		if not args.run_dir:
			parser.error("--format bout needs --run-dir")
		main_bout(args)


if __name__ == "__main__":
	main()
