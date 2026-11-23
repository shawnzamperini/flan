"""
Writes the synthetic W on W yield tables for reg_self_sputter: a constant
yield Y0 at every energy and angle, so every W atom or ion landing on the
wall sputters Y0 atoms whatever it hits with. Run it in the directory
reg_self_sputter runs in, before the y05 and y125 variants:

	python make_self_tables.py

It writes self_y05.h5 (Y0 = 0.5) and self_y125.h5 (Y0 = 1.25), in the
layout python/rustbca_yields.py writes and Flan reads. The lowest energy
is 1e-3 eV, so even slow neutrals sputter (the yield below a table's
lowest energy is 0).
"""

import argparse

import h5py
import numpy as np


TABLES = {"self_y05.h5": 0.5, "self_y125.h5": 1.25}


def write_table(path, y0):
	"""
	Writes a constant-yield W on W table.

	Inputs:
	 - path (str) : Output file.
	 - y0 (float) : The yield at every energy and angle.
	"""
	energy = np.array([1e-3, 1.0, 1e3, 1e6])
	angle = np.array([0.0, 30.0, 60.0, 89.0])
	with h5py.File(path, "w") as f:
		f["energy"] = energy
		f["angle"] = angle
		f["yield"] = np.full((len(energy), len(angle)), y0)
		f["projectile/Z"] = 74.0
		f["projectile/mass_amu"] = 183.84
		f["target/Z"] = 74.0
		f["target/mass_amu"] = 183.84
		f["target/Es"] = 8.79
		f.attrs["command"] = "regression/make_self_tables.py"
		f.attrs["description"] = "constant yield {} for reg_self_sputter" \
			.format(y0)


def main():
	parser = argparse.ArgumentParser(description=__doc__,
		formatter_class=argparse.RawDescriptionHelpFormatter)
	parser.parse_args()
	for path, y0 in TABLES.items():
		write_table(path, y0)
		print("Wrote {} (Y0 = {})".format(path, y0))


if __name__ == "__main__":
	main()
