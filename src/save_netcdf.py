# Save results into a netcdf file.
from netCDF4 import Dataset
from datetime import datetime


def save_netcdf(input_opts, gkyl_bkg, imp_dict):
	"""
	
	"""
	
	# Create NetCDF4 file with same name as input file.
	nc_fname = input_opts["case_name"] + ".nc"
	rootgrp = Dataset(nc_fname, "w")
	
	# Create groups for the input options, the Gkyell background and
	# the impurity transport results.
	inputgrp = rootgrp.createGroup("input_opts")
	gkylgrp = rootgrp.createGroup("gkyl_bkg")
	impgrp = rootgrp.createGroup("imp_results")
	
	# Give the groups some attributes. 
	time_str = datetime.now().strftime("%m/%d/%Y, %H:%M:%S")
	rootgrp.description = "Flan simulation results"
	rootgrp.history = "Run date: {}".format(time_str)
	inputgrp.description = "Input options for simulation"
	gkylgrp.description = "Background plasma from Gkeyll"
	impgrp.description = "Impurity transport results"

	# Create dimensions that will be assigned to the data.
	gkyl_xdim = rootgrp.createDimension("gkyl_xdim", len(gkyl_bkg["x"]))
	gkyl_ydim = rootgrp.createDimension("gkyl_ydim", len(gkyl_bkg["y"]))
	gkyl_zdim = rootgrp.createDimension("gkyl_zdim", len(gkyl_bkg["z"]))
	
	# Create the variables for the Gkeyll background.
	gkyl_x = gkylgrp.createVariable("gkyl_x", "f4", ("gkyl_xdim",))
	gkyl_y = gkylgrp.createVariable("gkyl_y", "f4", ("gkyl_ydim",))
	gkyl_z = gkylgrp.createVariable("gkyl_z", "f4", ("gkyl_zdim",))
	gkyl_ne = gkylgrp.createVariable("gkyl_ne", "f4", ("gkyl_xdim", 
		"gkyl_ydim", "gkyl_zdim"))
	
	# Add units to the variables.
	gkyl_x.units = "m"
	gkyl_y.units = "m"
	gkyl_z.units = "m"
	gkyl_ne.units = "m-3"
	
	# Put variables into NetCDF files.
	gkyl_x[:] = gkyl_bkg["x"]
	gkyl_y[:] = gkyl_bkg["y"]
	gkyl_z[:] = gkyl_bkg["z"]
	gkyl_ne[:] = gkyl_bkg["ne"]

	
