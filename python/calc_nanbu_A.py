# This script is a simple helper script to calculate the tabulated values for
# the variable A needed in the Nanbu collision model. A is solved from an 
# implicit equation, hence it is much computationally more efficient to 
# precalculate and then interpolate accordingly within Flan.
# I'm not gonna lie, AI wrote this code because I am running out of budget
# and time for this project :(

import numpy as np
from scipy.optimize import root_scalar


def equation(A, s):
	return np.cosh(A) / np.sinh(A) - 1/A - np.exp(-s)

def solve_for_A(s):
	
	# We apply approximations mentioned in Nanbu because outside this range
	# we actually struggle to converge anyways.
	if s < 0.01:
		return 1 / s
	elif s > 5:
		return 3 * np.exp(-s)
	else:
		result = root_scalar(equation, args=(s,), bracket=[-0.001, 300], 
			method='brentq')
		return result.root if result.converged else None

# Example usage
#s = 0.5
#A = solve_for_A(s)
#print(f"Solved A for s={s}: {A}")

# Print out to prove correct values. These are the value from Nanbu Table I. 
s_test = [0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.1, 0.2, 0.3,
	0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 2.0, 3.0, 4.0]
a_test = [100.5, 50.50, 33.84, 25.50, 20.50, 17.17, 14.79, 13.01, 11.62, 10.51, 
	5.516, 3.845, 2.987, 2.448, 2.067, 1.779, 1.551, 1.363, 1.207, 0.4105, 
	0.1496, 0.05496]
print("s     A (Nanbu)  A (calc)    Percent Error")
for i in range(len(s_test)):
	ans = solve_for_A(s_test[i])
	print("{:4.2f}  {:9.5f}  {:9.5f}  {:5.2f}%".format(s_test[i], a_test[i], ans,
		(ans - a_test[i]) / a_test[i] * 100)) 

# Tabulate A over 0.009 <= s <= 5.001. Rounding makes the s values exact
# (np.arange alone gives 0.009999... for 0.01, which falls below the s < 0.01
# cutoff in solve_for_A). Flan interpolates in this table assuming the s
# values are evenly spaced, so keep them that way.
s_values = np.round(np.arange(0.009, 5.002, 0.001), 3)
a_values = [solve_for_A(s) for s in s_values]

vals_per_line = 3

# Doxygen @file block written at the top of each header. The values are only
# written by this script, so describe them here instead of by hand.
def file_block(name, extra):
	return (
		"/**\n"
		f"* @file  {name}\n"
		"*\n"
		"* @brief Tabulated values of A(s) for the Nanbu collision model.\n"
		"*\n"
		"* Generated automatically by python/calc_nanbu_A.py, don't edit by\n"
		"* hand. A is found from coth(A) - 1/A = exp(-s) for evenly spaced s\n"
		"* from 0.009 to 5.001 in steps of 0.001. Outside 0.01 <= s <= 5,\n"
		"* Nanbu's approximations A = 1/s and A = 3 exp(-s) are used\n"
		"* instead.\n"
		f"{extra}"
		"*/\n")

def write_values(f, values, indent):
	"""
	Writes values to f, vals_per_line per line, separated by commas.
	"""
	for i in range(0, len(values), vals_per_line):
		chunk = values[i:i+vals_per_line]
		line = ", ".join("{:.6e}".format(v) for v in chunk)
		end = "," if i + vals_per_line < len(values) else ""
		f.write(f"{indent}{line}{end}\n")

# CPU version, include/nanbu_s_a.h
with open("nanbu_s_a.h", "w") as f:

	# Header information
	f.write("#ifndef NANBU_S_A_H\n")
	f.write("#define NANBU_S_A_H\n\n")
	f.write(file_block("nanbu_s_a.h", "*\n"
		"* CPU version. The GPU version is cuda/nanbu_s_a.cuh.\n"))
	f.write("#include <array>\n")
	f.write("namespace Nanbu\n{\n")

	# Write out s values
	f.write("std::array<double, {:}> s {{\n".format(len(s_values)))
	write_values(f, s_values, "    ")
	f.write("};\n\n")

	# Write out A values
	f.write("std::array<double, {:}> A {{\n".format(len(a_values)))
	write_values(f, a_values, "    ")
	f.write("};\n\n")

	# Closing namespace bracket
	f.write("}\n\n")
	f.write("#endif")

# GPU version, include/cuda/nanbu_s_a.cuh. Same values, but as plain
# __device__ arrays since std::array can't be used in device code.
with open("nanbu_s_a.cuh", "w") as f:

	# Header information
	f.write("#pragma once\n\n")
	f.write(file_block("nanbu_s_a.cuh", "*\n"
		"* GPU version, as plain __device__ arrays since std::array can't be\n"
		"* used in device code. The CPU version is nanbu_s_a.h.\n"))
	f.write("namespace Nanbu\n{\n")
	f.write("\t__device__ constexpr size_t N = {:};\n\n".format(len(s_values)))

	# Write out s values
	f.write("\t__device__ const double s[N] = {\n")
	write_values(f, s_values, "\t\t")
	f.write("\t};\n\n")

	# Write out A values
	f.write("\t__device__ const double A[N] = {\n")
	write_values(f, a_values, "\t\t")
	f.write("\t};\n\n")

	# Closing namespace bracket
	f.write("}\n")

print("nanbu_s_a.h and nanbu_s_a.cuh have been created. Move nanbu_s_a.h to")
print("include/ and nanbu_s_a.cuh to include/cuda/ so they get compiled in")
print("like any other header file.")
