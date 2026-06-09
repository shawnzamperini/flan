"""
Read an EQDSK/gfile equilibrium file into a Python dictionary.
Mirrors the Fortran read logic including optional sections (kvtor, nmass, keecur, iplcout).
"""

import re
import numpy as np

# Matches a floating-point number in scientific notation, including a leading sign.
# Handles run-together values like "4.167022815e-01-3.764044606e+00" where a negative
# number immediately follows a positive one with no whitespace separator.
_FLOAT_RE = re.compile(r"[+-]?\d+\.\d+[eEdD][+-]\d+|[+-]?\d+\.?\d*")


def _parse_floats_from_line(line):
	"""
	Tokenise one line of gfile float data.

	Handles both space-separated and run-together values such as:
		4.167022815e-01-3.764044606e+00
	where the minus sign of a negative number is the only delimiter.
	"""
	# Insert a space before any sign character that immediately follows an
	# exponent digit (i.e. the end of one number), then split normally.
	# Pattern: digit or sign, then look-ahead for [+-] not preceded by [eEdD].
	line = re.sub(r"(?<=[0-9])(?=[+-](?![eEdD]))", " ", line)
	return [float(x) for x in line.split()]


def _read_floats(f, n):
	"""Read exactly n floats from the file, spanning as many lines as needed."""
	vals = []
	while len(vals) < n:
		line = f.readline()
		if line == "":
			raise EOFError(f"Unexpected end of file while reading {n} floats")
		chunk = _parse_floats_from_line(line)
		vals.extend(chunk)
	if len(vals) != n:
		raise ValueError(f"Expected {n} floats, got {len(vals)}")
	return vals


def _read_ints_fmt2022(line):
	"""Format 2022: 2i5"""
	return int(line[0:5]), int(line[5:10])


def _read_fmt2024(line):
	"""Format 2024: i5, e16.9, i5"""
	kvtor  = int(line[0:5])
	rvtor  = float(line[5:21])
	nmass  = int(line[21:26])
	return kvtor, rvtor, nmass


def _read_fmt2026(line):
	"""Format 2026: i5"""
	return int(line[0:5])


def _read_fmt3000(line):
	"""Format 3000: 4i5"""
	return (int(line[0:5]), int(line[5:10]),
			int(line[10:15]), int(line[15:20]))


def _read_fmt3003(line):
	"""Format 3003: 2i5, i6, i5"""
	return (int(line[0:5]), int(line[5:10]),
			int(line[10:16]), int(line[16:21]))


def read_gfile(path, iplcout=0, nfcoil=0, nesum=0, ishot=0, itime=0):
	"""
	Read an EQDSK / gfile into a dictionary.

	Parameters
	----------
	path	 : str	– path to the gfile
	iplcout  : int	– controls optional PLC output section (default 0)
	nfcoil	 : int	– number of field coils (needed when iplcout == 1)
	nesum	 : int	– number of external currents (needed when iplcout == 1)
	ishot	 : int	– shot number (needed when iplcout == 1)
	itime	 : int	– time (needed when iplcout == 1)

	Returns
	-------
	dict with all quantities read from the file.
	"""
	g = {}

	with open(path, "r") as f:

		# ── Header line: format 2000  (6a8, 3i4) ────────────────────────────
		# Real files often have variable-width text; the last three tokens are
		# always idum, nw, nh as integers regardless of spacing.
		header = f.readline().rstrip("\n")
		tokens = header.split()
		idum, nw, nh = int(tokens[-3]), int(tokens[-2]), int(tokens[-1])
		text_part = header[:header.rfind(tokens[-3])].rstrip()
		case = [(text_part[i*8:(i+1)*8] if (i+1)*8 <= len(text_part)
				 else text_part[i*8:].ljust(8)) for i in range(6)]

		g["case"]  = case
		g["idum"]  = idum
		g["nw"]    = nw
		g["nh"]    = nh

		# ── Scalar quantities (format 2020: 5e16.9) ──────────────────────────
		v = _read_floats(f, 5)
		g["rdim"], g["zdim"], g["rcentr"], g["rleft"], g["zmid"] = v

		v = _read_floats(f, 5)
		g["rmaxis"], g["zmaxis"], g["simag"], g["sibry"], g["bcentr"] = v

		v = _read_floats(f, 5)
		g["current"] = v[0]
		g["simag"]	 = v[1]   # repeated
		g["rmaxis"]  = v[3]   # repeated (xdum at [2] and [4] discarded)

		v = _read_floats(f, 5)
		g["zmaxis"] = v[0]	  # repeated (xdum at [1],[3],[4] discarded)
		g["sibry"]	= v[2]	  # repeated

		# ── 1-D arrays of length nw ──────────────────────────────────────────
		g["fpol"]	= np.array(_read_floats(f, nw))
		g["pres"]	= np.array(_read_floats(f, nw))
		g["ffprim"] = np.array(_read_floats(f, nw))
		g["pprime"] = np.array(_read_floats(f, nw))

		# ── 2-D psi array: shape (nw, nh) ────────────────────────────────────
		# Fortran order: ((psirz(i,j),i=1,nw),j=1,nh)  →  C order (nh, nw) then transpose
		psirz_flat = _read_floats(f, nw * nh)
		g["psirz"] = np.array(psirz_flat).reshape(nh, nw).T   # shape (nw, nh)

		g["qpsi"] = np.array(_read_floats(f, nw))

		# ── Boundary & limiter counts (format 2022: 2i5) ─────────────────────
		line		  = f.readline()
		nbbbs, limitr = _read_ints_fmt2022(line)
		g["nbbbs"]	= nbbbs
		g["limitr"] = limitr

		# ── Boundary & limiter coordinates ───────────────────────────────────
		rb_zb		= _read_floats(f, 2 * nbbbs)
		g["rbbbs"]	= np.array(rb_zb[0::2])
		g["zbbbs"]	= np.array(rb_zb[1::2])

		rl_zl		= _read_floats(f, 2 * limitr)
		g["rlim"]	= np.array(rl_zl[0::2])
		g["zlim"]	= np.array(rl_zl[1::2])

		# ── Optional rotation / mass sections (format 2024: i5 e16.9 i5) ─────
		# All remaining sections are optional; many files end after the limiter.
		g["kvtor"] = 0;  g["rvtor"] = 0.0;	g["nmass"] = 0
		g["pressw"] = np.zeros(nw);  g["pwprim"] = np.zeros(nw)
		g["dmion"]	= np.zeros(nw)
		g["rhovn"]	= np.zeros(nw)
		g["keecur"] = 0
		g["epoten"] = np.zeros(nw)

		line = f.readline()
		if not line.strip():		  # EOF or blank → no optional sections
			pass
		else:
			kvtor, rvtor, nmass = _read_fmt2024(line)
			g["kvtor"] = kvtor
			g["rvtor"] = rvtor
			g["nmass"] = nmass

			if kvtor > 0:
				g["pressw"] = np.array(_read_floats(f, nw))
				g["pwprim"] = np.array(_read_floats(f, nw))

			if nmass > 0:
				g["dmion"] = np.array(_read_floats(f, nw))

			g["rhovn"] = np.array(_read_floats(f, nw))

			# ── keecur section (format 2026: i5) ─────────────────────────────
			line = f.readline()
			if line.strip():
				keecur = _read_fmt2026(line)
				g["keecur"] = keecur
				if keecur > 0:
					g["epoten"] = np.array(_read_floats(f, nw))

		# ── iplcout optional section ──────────────────────────────────────────
		g["iplcout"] = iplcout

		if iplcout > 0:
			if iplcout == 1:
				line = f.readline()
				if ishot <= 99999:
					nw2, nh2, ishot2, itime2 = _read_fmt3000(line)
				else:
					nw2, nh2, ishot2, itime2 = _read_fmt3003(line)
				g["ishot"] = ishot2
				g["itime"] = itime2

				v = _read_floats(f, 4)
				g["rgrid"] = np.zeros(nw);	g["rgrid"][0]  = v[0]; g["rgrid"][-1]  = v[1]
				g["zgrid"] = np.zeros(nh);	g["zgrid"][0]  = v[2]; g["zgrid"][-1]  = v[3]

				g["brsp"]	= np.array(_read_floats(f, nfcoil))
				g["ecurrt"] = np.array(_read_floats(f, nesum))
				g["pcurrt"] = np.array(_read_floats(f, nw * nh))

			elif iplcout == 2:
				pcurrz_flat = _read_floats(f, nw * nh)
				g["pcurrz"] = np.array(pcurrz_flat).reshape(nh, nw).T

				g["cjor"]	= np.array(_read_floats(f, nw))
				g["r1surf"] = np.array(_read_floats(f, nw))
				g["r2surf"] = np.array(_read_floats(f, nw))
				g["volp"]	= np.array(_read_floats(f, nw))
				g["bpolss"] = np.array(_read_floats(f, nw))

		# These are not included and must be created ourselves
		g["rgrid"] = np.linspace(g["rleft"], g["rleft"] + g["rdim"], g["nw"])
		g["zgrid"] = np.linspace(g["zmid"] - g["zdim"]/2, 
			g["zmid"] + g["zdim"]/2, g["nh"])

		# ── Uniform psi grid for all 1-D profile arrays ──────────────────────
		# fpol, pres, ffprim, pprime, and qpsi are all sampled on this grid.
		g["psi_grid"] = np.linspace(g["simag"], g["sibry"], nw)

	return g


# ── Example usage ─────────────────────────────────────────────────────────────
if __name__ == "__main__":
	import sys

	path = sys.argv[1] if len(sys.argv) > 1 else "g000000.00000"
	data = read_gfile(path)

	print(f"nw={data['nw']}, nh={data['nh']}")
	print(f"rmaxis={data['rmaxis']:.6f}, zmaxis={data['zmaxis']:.6f}")
	print(f"bcentr={data['bcentr']:.6f}, current={data['current']:.6e}")
	print(f"psirz shape: {data['psirz'].shape}")
	print(f"nbbbs={data['nbbbs']}, limitr={data['limitr']}")
	print("Keys:", sorted(data.keys()))

