# Script to calculate the electric field by calculating the gradient of the
# potential on an irregular grid. This takes advantage of the scipy.spatial
# library, hence why this python script exists because I am NOT going to try
# and implement that in C++. 
import numpy as np
from scipy.spatial import cKDTree
from tqdm import tqdm
from multiprocessing import Pool, Value
import argparse



# Global counter to keep track of progress across processes
counter = None

def load_XYZ():
	"""
	Load Cartesian X,Y,Z coordinates

	Returns
	-------
	XYZ : 3D numpy array (x,y,z)
		X,Y,Z coordinates as a 3D numpy array (x,y,z).
	"""

	# Name of file, shouldn't change
	XYZ_file = "bkg_from_pgkyl_cell_XYZ.csv"

	# Pretty simple. First 7 lines are comments, 8th line is an integer 
	# of how many entries follow (which isn't actually needed). The data is
	# printed flattened using C-style indexing.
	return np.loadtxt(XYZ_file, skiprows=8)


def load_bkg_file(fname):
	"""
	Returns data form a background file as a 4D numpy array [t,x,y,z]

	Parameters
	----------
	fname : str
		Filename of data that can be represented as 4D. E.g., 
		bkg_from_pgkyl_potential.csv

	Returns
	-------
	bkg_data : 4D np.array
		The background data in a numpy array
	"""

	# Line with data dimensions, also shouldn't change
	dim_line_num = 6

	# Need to load the shape of the data so we know how many frames there are.
	# Just read dim_line_num and break.
	with open(fname, "r") as f:
		for line_num, line in enumerate(f, start=1):
			if line_num == dim_line_num:
				
				# 4 int dimensions to load
				dims = np.array(line.split(), dtype=int)
				break

	# Read potential, reshape and return
	bkg_data = np.loadtxt(fname, skiprows=6).reshape(dims)
	return bkg_data


#def calc_gradients_3d(XYZ, pot, bmag):
def calc_gradients_3d(XYZ_pot_bmag_ntimes):
	"""
	Calculate potential and magnetic field gradients

	Parameters
	----------
	XYZ : 3D numpy array
		Contains the Cartesian X,Y,Z values at a given computational x,y,z

	pot : 4D numpy array (t,x,y,z)
		The plasma potential

	bmag : 4D numpy array (t,x,y,z)
		The magnetic field magnitude

	Returns
	-------
	eX, eY, eZ, gradB_X, gradB_Y, gradB_Z : 4D numpy arrays of the gradients
		Returns in this order as tuple
	"""
	# Want to use the global counter that is setup in shared memory
	global counter

	# Unpack arguments
	XYZ, pot, bmag, ntimes = XYZ_pot_bmag_ntimes

	# No t dimension, hence the [1:] here. Also append 3 to preserve coordinate
	# groupings.
	XYZ_shape = list(pot.shape[1:])
	XYZ_shape.append(3)
	XYZ_shaped = XYZ.reshape(XYZ_shape)

	def calculate_gradient(i, j, k, values):

		#print("Point ({},{},{})".format(i,j,k))
		#print(XYZ_shaped[i,j,k])
	
		# The 6 neighboring cells and a mask that starts all True. Edge cases
		# may assign False to mask elements.
		neighbors = [(i+1,j,k), (i,j+1,k), (i,j,k+1), (i-1,j,k), 
			(i,j-1,k), (i,j,k-1)]
		mask = [True, True, True, True, True, True]

		# Edge or grid cases to remove non-existant neighbors by assigning the 
		# mask to False
		if (i == 0): mask[3] = False 
		if (i == values.shape[0]-1): mask[0] = False
		if (j == 0): mask[4] = False 
		if (j == values.shape[1]-1): mask[1] = False
		if (k == 0): mask[5] = False 
		if (k == values.shape[2]-1): mask[2] = False

		# Remove any non-existing neighbors
		neighbors = [neighbors[m] for m in range(len(neighbors)) if mask[m]]

		# Skip the current point itself (indices[0] is the same point)
		neighbor_points = [XYZ_shaped[idx] for idx in neighbors]
		neighbor_values = [values[idx] for idx in neighbors]
		#print(neighbor_points)
		#print(neighbor_values)

		# Build the system of equations to solve (least squares)
		A = neighbor_points - XYZ_shaped[i,j,k]  # Differences in positions
		b = neighbor_values - values[i,j,k]  # Differences in scalar values

		# Solve the least squares problem A * grad = b
		grad, res, _, _ = np.linalg.lstsq(A, b, rcond=None)
		#print("grad & res")
		#print(grad)
		#print(res)

		return grad

	# Empty 4D (t,x,y,z) arrays for each electric field component
	elec_X = np.zeros(pot.shape)
	elec_Y = np.zeros(pot.shape)
	elec_Z = np.zeros(pot.shape)
	gradB_X = np.zeros(bmag.shape)
	gradB_Y = np.zeros(bmag.shape)
	gradB_Z = np.zeros(bmag.shape)

	# Go through one frame at a time
	# To-do: Parallelize this loop, probably with numba
	#print("Calculating potential and magnetic field gradients for each frame...")
	#for t in tqdm(range(pot.shape[0])):
	for t in range(pot.shape[0]):
		for i in range(pot.shape[1]):
			for j in range(pot.shape[2]):
				for k in range(pot.shape[3]):

					# Electric field
					grad = calculate_gradient(i, j, k, pot[t])	
					elec_X[t,i,j,k] = -grad[0]
					elec_Y[t,i,j,k] = -grad[1]
					elec_Z[t,i,j,k] = -grad[2]

					# Magnetic field gradient
					grad = calculate_gradient(i, j, k, bmag[t])	
					gradB_X[t,i,j,k] = grad[0]
					gradB_Y[t,i,j,k] = grad[1]
					gradB_Z[t,i,j,k] = grad[2]

		# Update shared counter, print progress update
		with counter.get_lock():
			counter.value += 1
			perc = int(counter.value / ntimes * 100)
			print("Frames: {}/{} ({:d}%)".format(counter.value, ntimes, perc), 
				end="\r", flush=True)

	# Return each component of electric field
	return elec_X, elec_Y, elec_Z, gradB_X, gradB_Y, gradB_Z

def write_field(field_XYZ, fname_base):
	"""
	Writes each component of the field to a file, e.g.,
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
		with open("{}_{}.csv".format(fname_base, comp), "w") as f:
			f.write(header)
			num_vals = "{:d} {:d} {:d} {:d}\n".format(field_XYZ[i].shape[0], 
					field_XYZ[i].shape[1], field_XYZ[i].shape[2], 
					field_XYZ[i].shape[3]) 
			f.write(num_vals)
			for j in range(0, len(field_XYZ[i])):
				np.savetxt(f, field_XYZ[i][j].flatten())

def init_counter(shared_counter):
	"""
	Function to be called at the beginning of each process to access
	shared_counter, which exists in shared memory.
	"""
	global counter
	counter = shared_counter

def main(num_processes=1):
	"""
	Controlling function for calculating electric field and magnetic field 
	gradients
	"""

	# Files needed: cell_center_XYZ.csv, bkg_from_pgkyl_potential.csv and
	# bkg_from_pgkyl_magnetic_magnitude.csv.
	XYZ = load_XYZ()
	pot = load_bkg_file("bkg_from_pgkyl_potential.csv")
	bmag = load_bkg_file("bkg_from_pgkyl_magnetic_magnitude.csv")

	# Make sure they are the same length. [0] here is because we have the
	# potential data for every frame, so just need any random frame to check
	# against.
	if (len(XYZ) != len(pot[0].flatten())):
		print("Error! The number of X,Y,Z coordinates does not match the ")
		print("number of plasma potential values!")
		print("  len(XYZ)={}  len(pot)={}".format(len(XYZ),	
			len(pot[0].flatten())))
		return None

	# Define chunk_size and break up into chunks (of frames) of that size. 
	# If num_processes > number of frames, then just reassign num_processes
	# to the number of frames.
	if (num_processes > pot.shape[0]): num_processes = pot.shape[0]
	chunk_size = pot.shape[0] // num_processes

	# Remainder if num_processes does not evenly divide the number of frames.
	remainder = pot.shape[0] % num_processes

	# Assemble chunks
	pot_chunks = []
	bmag_chunks = []
	for i in range(num_processes):
		start = i * chunk_size
		end = start + chunk_size

		# Last chunk absorbs remainder
		if i == num_processes -1:
			end += remainder

		pot_chunks.append(pot[start:end])
		bmag_chunks.append(bmag[start:end])

	#pot_chunks = [pot[i:i + chunk_size] for i in	
	#	range(0, pot.shape[0], chunk_size)]
	#bmag_chunks = [bmag[i:i + chunk_size] for i in	
	#	range(0, bmag.shape[0], chunk_size)]
	
	# If num_processes does not evenly go into the number of times, then we
	# need to assign the leftover frame to the last chunk.

	# Zip together into args to be passed below
	XYZs = [XYZ] * num_processes
	ntimes = [pot.shape[0]] * num_processes
	args = list(zip(XYZs, pot_chunks, bmag_chunks, ntimes))

	# Create process pool and get each process a chunk of frames
	shared_counter = Value("i", 0)
	#print("Calculating gradients...")
	with Pool(processes=num_processes, initializer=init_counter, 
		initargs=(shared_counter,)) as pool:

		# Calculate electric field and gradB components
		results = pool.imap(calc_gradients_3d, args)

		# Put each gradient result into separate lists
		eX_list, eY_list, eZ_list, gbX_list, gbY_list, gbZ_list  \
			= zip(*results)

	# Clean up terminal from progress printing
	print(" " * 50, end="\r")
	
	# Concatenate along time axis
	eX = np.concatenate(eX_list, axis=0)
	eY = np.concatenate(eY_list, axis=0)
	eZ = np.concatenate(eZ_list, axis=0)
	gbX = np.concatenate(gbX_list, axis=0)
	gbY = np.concatenate(gbY_list, axis=0)
	gbZ = np.concatenate(gbZ_list, axis=0)

	# Package up together in same format as if not multithreaded
	grads = (eX, eY, eZ, gbX, gbY, gbZ)
			
	# Separate out to be written
	elec_field_XYZ = grads[0:3]
	gradB_XYZ = grads[3:]

	# Write out the three X,Y,Z files each
	write_field(elec_field_XYZ, "bkg_from_pgkyl_elec_field")
	write_field(gradB_XYZ, "bkg_from_pgkyl_gradb")

if __name__ == "__main__":

	# Command line argument to access number of processes to spawn
	desc = "Calculate electric field and magentic field gradient."
	parser = argparse.ArgumentParser(description=desc)
	parser.add_argument("-n", "--num-processes", type=int, default=1,
		help="number of processes to run with")
	
	# Parse, defaults to one process if nothing provided
	args = parser.parse_args()

	# Run
	main(args.num_processes)

