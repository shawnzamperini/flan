# Script to calculate the electric field by calculating the gradient of the
# potential on an irregular grid. This takes advantage of the scipy.spatial
# library, hence why this python script exists because I am NOT going to try
# and implement that in C++. 
import numpy as np
from scipy.spatial import cKDTree
from tqdm import tqdm


def load_XYZ():
	"""
	Returns the X,Y,Z coordinates as a 3D numpy array [x,y,z].
	"""

	# Name of file, shouldn't change
	XYZ_file = "cell_center_XYZ.csv"
	"""
	# Line with data dimensions, also shouldn't change
	dim_line_num = 8

	# Need to load the shape of the data so we know how many frames there are.
	# Just read dim_line_num and break.
	with open(XYZ_file, "r") as f:
		for line_num, line in enumerate(f, start=1):
			if line_num == dim_line_num:
				
				# 4 int dimensions to load
				dims = np.array(line.split(), dtype=int)
				break

	"""
	# Pretty simple. First 7 lines are comments, 8th line is an integer 
	# of how many entries follow (which isn't actually needed). The data is
	# printed flattened using C-style indexing.
	return np.loadtxt(XYZ_file, skiprows=8)

def load_potential():
	"""
	Returns the potential as a 4D numpy array [t,x,y,z]
	"""

	# Name of file, shouldn't change
	pot_file = "bkg_from_pgkyl_potential.csv"

	# Line with data dimensions, also shouldn't change
	dim_line_num = 6

	# Need to load the shape of the data so we know how many frames there are.
	# Just read dim_line_num and break.
	with open(pot_file, "r") as f:
		for line_num, line in enumerate(f, start=1):
			if line_num == dim_line_num:
				
				# 4 int dimensions to load
				dims = np.array(line.split(), dtype=int)
				break

	# Read potential, reshape and return
	pot = np.loadtxt(pot_file, skiprows=6).reshape(dims)
	return pot

def calc_elec_field(XYZ, pot):
	"""
	Uses a k-nearest neighbor tree to calculate the electric field from the
	gradient of the potential on an irregular grid. Returns each X,Y,Z 
	component of the electric field, [elec_field_X, elec_field_Y, elec_field_Z],
	where each is a 4D numpy array of [t, x, y, z]. 
	"""
	
	# Array to hold electric field values
	elec_X = np.zeros(pot.shape)
	elec_Y = np.zeros(pot.shape)
	elec_Z = np.zeros(pot.shape)

	# This is a nearest neighbor tree, it is used to find the nearest
	# coordinates to calculate the gradient from
	tree = cKDTree(XYZ)

	# Function to calculate gradient at a specific point. Thanks Copilot!
	def calculate_gradient(idx, values):
		point = XYZ[idx]

		# Query the nearest neighbors. We choose the 6 nearest neighbors in
		# the spirit of getting the two nearest in each coordinate direction.
		# Note: getting the 6 nearest is k=7.
		distances, indices = tree.query(point, k=7)
		
		# Skip the current point itself (indices[0] is the same point)
		neighbors = indices[1:]
		neighbor_points = XYZ[neighbors]
		neighbor_values = values[neighbors]

		# Build the system of equations to solve (least squares)
		A = neighbor_points - point  # Differences in positions
		b = neighbor_values - values[idx]  # Differences in scalar values

		# Solve the least squares problem A * grad = b
		grad, _, _, _ = np.linalg.lstsq(A, b, rcond=None)
		return grad
	
	# Go through one frame at a time
	print("Calculating potential gradient for each frame...")
	for t in tqdm(range(0, pot.shape[0])):

		# The points are already in a 1D array of 3-values (X,Y,Z), so flatten
		# the potential values for this frame
		pot_vals = pot[t].flatten()
		elec_vals_X = np.zeros(pot_vals.shape)
		elec_vals_Y = np.zeros(pot_vals.shape)
		elec_vals_Z = np.zeros(pot_vals.shape)

		# Need to loop through every point and calculate the gradient there
		for i in range(0, len(pot_vals)):
			grad = calculate_gradient(i, pot_vals)
			elec_vals_X[i] = -grad[0]
			elec_vals_Y[i] = -grad[1]
			elec_vals_Z[i] = -grad[2]

		# Reshape calculate electric field to the spatial dimensions of the
		# grid, and then store for this frame.
		elec_X[t] = elec_vals_X.reshape(pot[t].shape)
		elec_Y[t] = elec_vals_Y.reshape(pot[t].shape)
		elec_Z[t] = elec_vals_Z.reshape(pot[t].shape)
	
	# Return each component of electric field
	return elec_X, elec_Y, elec_Z

def write_elec_field(elec_XYZ):
	"""
	Writes each component of the electric field to a file, e.g.,
	bkg_from_pgkyl_elec_field_X.csv.
	"""

	# Add an informative header just for clarity's sake
	header = (
			 "# The values for the specified type of data requested from\n"
			 "# Gkeyll. The first row are integers: the number of frames,\n"
			 "# and then the shape of each array (x, y, z dimensions).\n"
			 "# Then each array for each frame is printed out flattened\n"
			 "# using C-style indexing.\n"
			 ) 
	
	# Write out. i selects the X, Y or Z component, so we need a file for each
	# component
	for i, comp in enumerate(["X", "Y", "Z"]):
		with open("bkg_from_pgkyl_elec_field_{}.csv".format(comp), "w") as f:
			f.write(header)
			num_vals = "{:d} {:d} {:d} {:d}\n".format(elec_XYZ[i].shape[0], 
					elec_XYZ[i].shape[1], elec_XYZ[i].shape[2], 
					elec_XYZ[i].shape[3]) 
			f.write(num_vals)
			for j in range(0, len(elec_XYZ[i])):
				np.savetxt(f, elec_XYZ[i][j].flatten())

def main():
	"""
	Controlling function for calculating electric field
	"""

	# Files needed: cell_center_XYZ.csv and bkg_from_pgkyl_potential.csv
	XYZ = load_XYZ()
	pot = load_potential()

	# Make sure they are the same length. [0] here is because we have the
	# potential data for every frame, so just need any random frame to check
	# against.
	if (len(XYZ) != len(pot[0].flatten())):
		print("Error! The number of X,Y,Z coordinates does not match the ")
		print("number of plasma potential values!")
		print("  len(XYZ)={}  len(pot)={}".format(len(XYZ),	
			len(pot[0].flatten())))
		return None

	# Calculate electric field components from potential gradient
	elec_field_XYZ = calc_elec_field(XYZ, pot)

	# Write out to elec_field.csv
	write_elec_field(elec_field_XYZ)

if __name__ == "__main__":
	main()

