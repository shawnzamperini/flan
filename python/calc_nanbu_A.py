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

def solve_for_A(s, guess=1.0):
	
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

# Loop through a range of values, creating the file tabulated values nanbu_a.h
s_values = np.arange(0.009, 5.002, 0.001)
#s_values = np.arange(0.01, 5, 0.1)
a_values = [solve_for_A(s) for s in s_values]

vals_per_line = 3
with open("nanbu_s_a.h", "w") as f:

	# Header information
	f.write("#ifndef NANBU_S_A_H\n")
	f.write("#define NANBU_S_A_H\n\n")
	f.write("// This file generated automatically by python/calc_nanbu_A.py\n")
	f.write("#include <array>\n")
	f.write("namespace Nanbu\n{\n")

	# Write out s values
	f.write("std::array<double, {:}> s {{\n".format(len(s_values)))
	for i in range(0, len(s_values), vals_per_line):
		chunk = s_values[i:i+vals_per_line]
		line = ", ".join("{:.3f}".format(v) for v in chunk)
		end = "," if i + 3 < len(s_values) else ""
		f.write(f"    {line}{end}\n")
	f.write("};\n\n")

	# Write out A values
	f.write("std::array<double, {:}> A {{\n".format(len(a_values)))
	for i in range(0, len(a_values), vals_per_line):
		chunk = a_values[i:i+vals_per_line]
		line = ", ".join("{:.3f}".format(v) for v in chunk)
		end = "," if i + 3 < len(a_values) else ""
		f.write(f"    {line}{end}\n")
	f.write("};\n\n")

	# Closing namespace bracket
	f.write("}\n\n")
	f.write("#endif")

print("nanbu_s_a.h has been created. Move to include/ so it gets compiled in")
print("like any other header file.")
