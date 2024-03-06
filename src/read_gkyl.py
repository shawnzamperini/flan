# Read frames from a Gkeyll run. This file is subject to change, as the
# output of Gkeyll seems to still be in flux a little. 
import sys
import postgkyl as pg
import numpy as np
from tqdm import tqdm
import constants


def read_gkyl(input_opts):
	"""
	Top level function to interface with the different ways of loading
	Gkeyll data. 
	"""
	
	# Get the needed information from the input file.
	gkyl_dir = input_opts["gkyl_dir"]
	gkyl_name = input_opts["gkyl_name"]
	gkyl_format = input_opts["gkyl_format"]
	gkyl_fstart = input_opts["gkyl_fstart"]
	gkyl_fend = input_opts["gkyl_fend"]
	gkyl_elc_name = input_opts["gkyl_elc_name"]
	gkyl_ion_name = input_opts["gkyl_ion_name"]
	gkyl_bmag = input_opts["func_bmag"]
	
	# Select the type of Gkeyll data to load, which depends on the 
	# version of Gkeyll that was used. All will return gkyl_bkg of
	# which is a dictionary containing the background plasma from
	# Gkeyll that consists of: ne, te, ni, ti and elec. 
	if gkyl_format == "adios1":
		gkyl_bkg = read_adios1(gkyl_dir, gkyl_name, gkyl_fstart, 
			gkyl_fend, gkyl_elc_name, gkyl_ion_name, gkyl_bmag)
	elif gkyl_format() == "adios2":
		gkyl_bkg = read_adios2()
	elif gkyl_format == "gkylzero":
		gkyl_bkg = read_gkylzero()
	else:
		print("Error: Incorrect input for gkyl_format (). Valid " \
			"options are 'adios1', 'adios2' or 'gkylzero'."
			.format(gkyl_format))
		sys.exit()
	
	# Average all the frames together for steady-state background.
	gkyl_bkg["ne_avg"] = gkyl_bkg["ne"].mean(axis=0)
	gkyl_bkg["te_avg"] = gkyl_bkg["te"].mean(axis=0)
	gkyl_bkg["ni_avg"] = gkyl_bkg["ni"].mean(axis=0)
	gkyl_bkg["ti_avg"] = gkyl_bkg["ti"].mean(axis=0)
	
	
	return gkyl_bkg

def read_adios1(gkyl_dir, gkyl_name, fstart, fend, gkyl_elc_name, 
	gkyl_ion_name, gkyl_bmag):
	"""
	Read Gkeyll data that was created with ADIOS1. This is easily
	identified by if the .bp files are standalone (ADIOS1) or put into
	their own directories (ADIOS2, see other function). 
	
	gkyl_dir (str): Path to the directory with the .bp file.
	gkyl_name (str): Name of the Gkyell run.
	fstart (int): First frame of the Gkeyll simulation to load.
	fend (int): Last frame of the Gkeyll simulation to load.
	"""
	
	def read_bp_file(gkyl_param, gkyl_idx, gkyl_species=None, 
		multiplier=1.0):
		"""
		Function to read in a .bp file. 
		
		gkyl_param (str): Name of the Gkyell parameter to load. This 
			generally is one of "M0" (density), "Temp" or "phi". 
		gkyl_idx (int): Index of the .bp file to load (e.g., 
			my_gkyl_run_electrons_M0_56.bp, where 56 is gkyl_idx). 
		gkyl_species (str): Name of the species (e.g., 
			my_gkyl_run_electrons_M0_56.bp, where electrons is species).
			Leave this as None if there is no species. E.g., use None
			to load the potential (my_gkyl_run_phi_56.bp).
		multiplier (float): Some data (e.g., Temp) needs to be
			multiplied by a scaler to bring it into normal units. 
		"""
		
		# Create path based on if species is included or not. 
		if gkyl_species is None:
			path = "{}/{}_{}_{}.bp".format(gkyl_dir, gkyl_name, 
				gkyl_param, gkyl_idx)
		else:
			path = "{}/{}_{}_{}_{}.bp".format(gkyl_dir, gkyl_name, 
				gkyl_species, gkyl_param, gkyl_idx)
			
		# Load GData object, pull out some simulation parameters. 
		gdata = pg.GData(path)
		time = float(gdata.ctx["time"])
		poly_order = gdata.ctx["polyOrder"]
		basis_type = gdata.ctx["basisType"]
		# charge = gdata.ctx["charge"]
		# mass = gdata.ctx["mass"]
		modal_flag = gdata.ctx["isModal"]
		
		# Interpolate the data if modal. I've not encountered what to
		# do if this isn't the case, so exit if so.
		if modal_flag:
			gdata_interp = pg.GInterpModal(gdata, poly_order, 'ms')
		else:
			print("Error: isModal is False. Not sure how to handle" \
				" this. Exiting...")
			sys.exit()
			
		# Extract the 3D data. I admit I don't fully know what these
		# steps are, I'm just copying the website example. 
		grid, data = gdata_interp.interpolate(0)
		data = data[:, :, :, 0] * multiplier
		
		return time, grid, data
	
	# Create a dictionary which will hold the Gkeyll background data.
	gkyl_bkg = {k:[] for k in ["time", "ne", "te", "ni", "ti", "elecx", 
		"elecy", "elecz"]}
	
	# Grab all the Gkeyll data for the desired range of frames and 
	# consolidate into numpy arrays.
	print("Loading data from {}...".format(gkyl_dir))
	for f in tqdm(range(fstart, fend)):
		
		# Electron density and temperature.
		t, grid, tmp_ne = read_bp_file("M0", f, gkyl_elc_name)
		t, grid, tmp_te = read_bp_file("Temp", f, gkyl_elc_name, 
			1/constants.ev)
		
		# Ion density and temperature.
		t, grid, tmp_ni = read_bp_file("M0", f, gkyl_ion_name)
		t, grid, tmp_ti = read_bp_file("Temp", f, gkyl_ion_name, 
			1/constants.ev)
		
		# Plasma potential. 
		t, grid, tmp_phi = read_bp_file("phi", f)
		
		# Parallel deuterium flow speed (Upar). 
		t, grid, tmp_upar = read_bp_file("Upar", f, gkyl_ion_name)
	
		# Save the x, y, z values of the grid. Only need to do this 
		# once. Messy code, but I'm just copy/pasting and thus, lazy.
		if f == fstart:
			x = grid[0]
			dx = np.diff(x)[0]
			x = x + dx / 2
			x = x[:-1]

			y = grid[1]
			dy = np.diff(y)[0]
			y = y + dy / 2
			y = y[:-1]

			z = grid[2]
			dz = np.diff(z)[0]
			z = z + dz / 2
			z = z[:-1]
			
			gkyl_bkg["x"] = x
			gkyl_bkg["y"] = y
			gkyl_bkg["z"] = z
		
		# Electric field.
		tmp_elec = np.gradient(-tmp_phi, x, y, z)
	
		# Store everything into our dictionary.
		gkyl_bkg["time"].append(t)
		gkyl_bkg["ne"].append(tmp_ne)
		gkyl_bkg["te"].append(tmp_te)
		gkyl_bkg["ni"].append(tmp_ni)
		gkyl_bkg["ti"].append(tmp_ti)
		gkyl_bkg["elecx"].append(tmp_elec[0])
		gkyl_bkg["elecy"].append(tmp_elec[1])
		gkyl_bkg["elecz"].append(tmp_elec[2])
	
	# Convert everything to numpy arrays. We do this here because
	# appending to a list in the previous section is faster than
	# appending to a numpy array, but ultimately we want to work with 
	# numpy arrays.
	gkyl_bkg = {k:np.array(v) for k, v in gkyl_bkg.items()}
		
	# Fill in an array using the magnetic field function. Right now
	# assume a magnetic field that only depends on the x value, which
	# is assumed to be it's R value (this is the simple helical 
	# approximation). This of course needs to be generalized. 
	b = np.zeros(gkyl_bkg["ne"].shape[1:])
	for i in range(0, len(x)):
		b[i, :, :] = gkyl_bmag(x[i])
	
	# Calculate the magnetic field gradient, then store into gkyl_bkg.
	grad_b = np.gradient(b, x, y, z)
	gkyl_bkg["b"] = b
	gkyl_bkg["grad_bx"] = grad_b[0]
	gkyl_bkg["grad_by"] = grad_b[1]
	gkyl_bkg["grad_bz"] = grad_b[2]
	
	return gkyl_bkg


def read_adios2():
	"""
	Read Gkeyll data that was created with ADIOS2. This is easily
	identified by if the .bp files are standalone (ADIOS1, see other 
	function) or put into their own directories (ADIOS2). 
	"""
	pass
	
def read_gkylzero():
	"""
	Read Gkeyll data that was made with gkylzero. These files end in 
	.gkyl, though this may be subject to change as gkylzero is 
	developed.
	"""
	pass
