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
	gkyl_num_interp = input_opts["gkyl_num_interp"]
	
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
	
	# Inteprolate extra frames between each frame if desired.
	if gkyl_num_interp is not None:
		gkyl_bkg = interp_frames(gkyl_bkg, gkyl_num_interp)

	# Average all the frames together for steady-state background.
	gkyl_bkg["ne_avg"] = gkyl_bkg["ne"].mean(axis=0)
	gkyl_bkg["te_avg"] = gkyl_bkg["te"].mean(axis=0)
	gkyl_bkg["ni_avg"] = gkyl_bkg["ni"].mean(axis=0)
	gkyl_bkg["ti_avg"] = gkyl_bkg["ti"].mean(axis=0)
	
	# Convert everything to floats. We use floats in impurity_following,
	# I'm sure this will bite me in the ass one day.
	for k, v in gkyl_bkg.items():
		gkyl_bkg[k] = gkyl_bkg[k].astype(np.float32)
	
	return gkyl_bkg

def read_adios1(gkyl_dir, gkyl_name, fstart, fend, gkyl_elc_name, 
	gkyl_ion_name, gkyl_bmag):
	"""
	Read Gkeyll data that was created with ADIOS1. This is easily
	identified by if the .bp files are standalone (ADIOS1) or put into
	their own directories (ADIOS2, see other function). 
	
	---Input---
	gkyl_dir (str): Path to the directory with the .bp file.
	gkyl_name (str): Name of the Gkyell run.
	fstart (int): First frame of the Gkeyll simulation to load.
	fend (int): Last frame of the Gkeyll simulation to load.
	gkyl_elc_name (str): The name used for the electron species.
	gkyl_ion_name (str): The name used for the (dueterium) ions.
	gkyl_bmag (func): This is a function that is input in the input file
		and then passed to this function. It contains the magnetic field
		as a function of the x, y, and z coordinates.
	"""
	
	def read_bp_file(gkyl_param, gkyl_idx, gkyl_species=None, 
		multiplier=1.0):
		"""
		Function to read in a .bp file. 
		
		---Input---
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
		"elecy", "elecz", "dtedx", "dtedy", "dtedz", "dtidx", "dtidy", 
		"dtidz", "viz"]}
	
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
		
		# Temperature gradients.
		tmp_gradte = np.gradient(tmp_te, x, y, z)
		tmp_gradti = np.gradient(tmp_ti, x, y, z)
	
		# Store everything into our dictionary.
		gkyl_bkg["time"].append(t)
		gkyl_bkg["ne"].append(tmp_ne)
		gkyl_bkg["te"].append(tmp_te)
		gkyl_bkg["ni"].append(tmp_ni)
		gkyl_bkg["ti"].append(tmp_ti)
		gkyl_bkg["elecx"].append(tmp_elec[0])
		gkyl_bkg["elecy"].append(tmp_elec[1])
		gkyl_bkg["elecz"].append(tmp_elec[2])
		gkyl_bkg["dtedx"].append(tmp_gradte[0])
		gkyl_bkg["dtedy"].append(tmp_gradte[1])
		gkyl_bkg["dtedz"].append(tmp_gradte[2])
		gkyl_bkg["dtidx"].append(tmp_gradti[0])
		gkyl_bkg["dtidy"].append(tmp_gradti[1])
		gkyl_bkg["dtidz"].append(tmp_gradti[2])
		gkyl_bkg["viz"].append(tmp_upar)
	
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
	gkyl_bkg["b"] = b.astype(np.single)
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

def interp_frames(gkyl_bkg, gkyl_num_interp):
	"""
	Interpolate a number of frames between each Gkeyll frame for
	each plasma quantity.
	
	---Input---
	gkyl_bkg (dict): The Gkeyll background returns from one of the read 
		functions above.
	gkyl_num_interp (int): The number of frames to add between each 
		Gkeyll frame.
		
	---Output---
	gkyl_bkg (dict): The Gkeyll background in the same format as before,
		just with the additional frames added in.
	
	"""

	# Extract how many frames we currently have, as well as the timestep
	# both before and after the interpolating operation.
	nframes = gkyl_bkg["ne"].shape[0]
	old_dt = gkyl_bkg["time"][1] - gkyl_bkg["time"][0]
	new_dt = old_dt / (gkyl_num_interp + 1)
	new_nframes = nframes + (nframes-1) * gkyl_num_interp
	
	# Extract all the plasma parameters into a list of arrays to loop 
	# through.
	print("Interpolating frames...")
	print("Original timestep: = {:.2e}".format(old_dt))
	print("New timestep = {:.2e}".format(new_dt))
	keys = ["ne", "te", "ni", "ti", "elecx", "elecy", "elecz", "dtedx", 
		"dtedy", "dtedz", "dtidx", "dtidy", "dtidz", "viz"]
	arrs = {k: gkyl_bkg[k] for k in keys}
	
	# Loop through one array at a time, tmp_arrs will hold the 
	# interpolated arrays, which we will overwrite the old arrays with 
	# at the end.
	tmp_arrs = {k: None for k in keys}
	for k in tqdm(arrs.keys()):
		
		# Create temporary array of the right size, filled with zeros.
		# Note that we are just increasing the first dimension here (the
		# number of frames), the others stay the same since we are just
		# interpolating in time. In theory one could interpolate in 
		# space with the same method, but if that's the case probably
		# best to rerun the Gkeyll simulation.
		# The new number of frames is:
		#   nframes + (nframes-1) * gkyl_num_interp
		# The below picture may help understand, where * = orginal frame
		# and - = interpolated frame. Say nframes = 5 and 
		# gkyl_num_interp = 3. There are then four gaps (nframes-1) to 
		# interpolate three frames into. The new number of frames is the
		# original number plus 5 + (5-1) * 3 = 17.
		#
		#    0   1   2   3   4
		#    *   *   *   *   *      I've 0-indexed these because
		#    *---*---*---*---*      of python
		#    0123456789...
		#
		# We can also infer the old-->new mapping for frame number here.
		#   0-->0
		#   1-->4
		#   2-->8
		#   3-->12
		#   4-->16
		# So the mapping is to just multiply by (gkyl_num_interp+1).
		tmp_arr = np.zeros((new_nframes, arrs[k].shape[1], 
			arrs[k].shape[2], arrs[k].shape[3]))
		for f in range(0, nframes-1):
			
			# Just use good ole point-slope formula for the interpolated 
			# frames.
			y0 = arrs[k][f]
			y1 = arrs[k][f + 1]
			m = (y1 - y0) / old_dt
			tmp_arr[f * (gkyl_num_interp+1)] = y0
			tmp_arr[(f+1) * (gkyl_num_interp+1)] = y1
			for g in range(1, gkyl_num_interp+1):
				tmp_arr[f * (gkyl_num_interp+1)+g] = m * new_dt * g + y0
		tmp_arrs[k] = tmp_arr
	
	# Reassign into gkyl_bkg.
	for k in keys:
		gkyl_bkg[k] = tmp_arrs[k] 
		
	# Don't forget to update the time array.
	old_tstart = gkyl_bkg["time"][0]
	old_tend = gkyl_bkg["time"][-1]
	gkyl_bkg["time"] = np.arange(old_tstart, old_tend+new_dt, new_dt)
	
	return gkyl_bkg
