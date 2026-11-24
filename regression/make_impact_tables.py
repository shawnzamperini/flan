"""
Writes the synthetic D on W yield tables for reg_sheath_impacts, flat in
angle so the sputtered source only depends on the D+ impact energies. Run
it in the directory reg_sheath_impacts runs in:

	python make_impact_tables.py

It writes, in the layout python/rustbca_yields.py writes and Flan reads,

 - impact_const.h5 : Y = 0.02 at every energy and angle, so the source
                     doesn't depend on the impacts at all.
 - impact_lin.h5   : Y = 1e-3 E / eV, so the yield averaged over the
                     impacts is 1e-3 times their mean energy. Flan
                     interpolates linearly in log10(E), so the table has
                     100 energies a decade (from 1e-3 to 1e5 eV) to keep
                     the interpolation error of a linear yield below 1e-4.
"""

import argparse

import h5py
import numpy as np


LIN_SLOPE = 1e-3    # yield per eV of impact_lin.h5


def write_table(path, energy, yld, description):
	"""
	Writes a D on W table that is flat in angle.

	Inputs:
	 - path (str)           : Output file.
	 - energy (array)       : Energies (eV).
	 - yld (array)          : Yield at each energy.
	 - description (str)    : Saved as an attribute.
	"""
	angle = np.array([0.0, 30.0, 60.0, 89.0])
	with h5py.File(path, "w") as f:
		f["energy"] = energy
		f["angle"] = angle
		f["yield"] = np.repeat(yld[:, None], len(angle), axis=1)
		f["projectile/Z"] = 1.0
		f["projectile/mass_amu"] = 2.014
		f["target/Z"] = 74.0
		f["target/mass_amu"] = 183.84
		f["target/Es"] = 8.79
		f.attrs["command"] = "regression/make_impact_tables.py"
		f.attrs["description"] = description


def main():
	parser = argparse.ArgumentParser(description=__doc__,
		formatter_class=argparse.RawDescriptionHelpFormatter)
	parser.parse_args()

	energy = np.array([1e-3, 1.0, 1e3, 1e6])
	write_table("impact_const.h5", energy, np.full(len(energy), 0.02),
		"constant yield 0.02 for reg_sheath_impacts")
	print("Wrote impact_const.h5 (Y = 0.02)")

	energy = np.logspace(-3.0, 5.0, 801)
	write_table("impact_lin.h5", energy, LIN_SLOPE * energy,
		"yield {} E / eV for reg_sheath_impacts".format(LIN_SLOPE))
	print("Wrote impact_lin.h5 (Y = {} E / eV)".format(LIN_SLOPE))


if __name__ == "__main__":
	main()
