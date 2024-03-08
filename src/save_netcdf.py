# Save results into a netcdf file.
from netCDF4 import Dataset
from datetime import datetime
import numpy as np


def save_netcdf(input_opts, gkyl_bkg, imp_dict, compression="zlib"):
	"""
	
	"""
	
	# Create NetCDF4 file with same name as input file.
	print("NetCDF: Creating file...")
	nc_fname = input_opts["case_name"] + ".nc"
	rootgrp = Dataset(nc_fname, "w")
	
	# Create groups for the input options, the Gkyell background and
	# the impurity transport results.
	inputgrp = rootgrp.createGroup("input_opts")
	gkylgrp = rootgrp.createGroup("gkyl_bkg")
	impgrp = rootgrp.createGroup("imp_results")
	
	# Give the groups some attributes. 
	time_str = datetime.now().strftime("%m/%d/%Y, %H:%M:%S")
	rootgrp.description = "Flan simulation results for case: {:}" \
		.format(input_opts["case_name"])
	rootgrp.history = "Run date: {}".format(time_str)
	inputgrp.description = "Input options for simulation"
	gkylgrp.description = "Background plasma from Gkeyll"
	impgrp.description = "Impurity transport results"

	# Create dimensions that will be assigned to the data.
	gkyl_xdim = rootgrp.createDimension("gkyl_xdim", len(gkyl_bkg["x"]))
	gkyl_ydim = rootgrp.createDimension("gkyl_ydim", len(gkyl_bkg["y"]))
	gkyl_zdim = rootgrp.createDimension("gkyl_zdim", len(gkyl_bkg["z"]))
	str_dim = rootgrp.createDimension("str_dim", 1)
	
	# In this case, the number of frames may not be the same at that 
	# specified in the input file (gkyl_fend-gkyl_fstart+1) if we 
	# interpolated frames. SO instead pull the number of frames from the
	# ne array.
	gkyl_nframes = gkyl_bkg["ne"].shape[0]
	gkyl_fdim = rootgrp.createDimension("gkyl_fdim", gkyl_nframes)
	
	# Create the variables for the Gkeyll background. We apply
	# compression as well as taking advantage of the significant_digits
	# keyword to cut down on storage.
	gkyl_x = gkylgrp.createVariable("gkyl_x", "f4", ("gkyl_xdim",), 
		compression=compression)
	gkyl_y = gkylgrp.createVariable("gkyl_y", "f4", ("gkyl_ydim",), 
		compression=compression)
	gkyl_z = gkylgrp.createVariable("gkyl_z", "f4", ("gkyl_zdim",), 
		compression=compression)
	gkyl_ne = gkylgrp.createVariable("gkyl_ne", "f4", ("gkyl_fdim", 
		"gkyl_xdim", "gkyl_ydim", "gkyl_zdim"), 
		compression=compression, significant_digits=6)
	gkyl_te = gkylgrp.createVariable("gkyl_te", "f4", ("gkyl_fdim", 
		"gkyl_xdim", "gkyl_ydim", "gkyl_zdim"), 
		compression=compression, significant_digits=6)
	gkyl_ni = gkylgrp.createVariable("gkyl_ni", "f4", ("gkyl_fdim", 
		"gkyl_xdim", "gkyl_ydim", "gkyl_zdim"), 
		compression=compression, significant_digits=6)
	gkyl_ti = gkylgrp.createVariable("gkyl_ti", "f4", ("gkyl_fdim", 
		"gkyl_xdim", "gkyl_ydim", "gkyl_zdim"), 
		compression=compression, significant_digits=6)
	gkyl_ex = gkylgrp.createVariable("gkyl_elecx", "f4", ("gkyl_fdim", 
		"gkyl_xdim", "gkyl_ydim", "gkyl_zdim"), 
		compression=compression, significant_digits=6)
	gkyl_ey = gkylgrp.createVariable("gkyl_elecy", "f4", ("gkyl_fdim", 
		"gkyl_xdim", "gkyl_ydim", "gkyl_zdim"), 
		compression=compression, significant_digits=6)
	gkyl_ez = gkylgrp.createVariable("gkyl_elecz", "f4", ("gkyl_fdim", 
		"gkyl_xdim", "gkyl_ydim", "gkyl_zdim"), 
		compression=compression, significant_digits=6)
	gkyl_b = gkylgrp.createVariable("gkyl_b", "f4", ( 
		"gkyl_xdim", "gkyl_ydim", "gkyl_zdim"), 
		compression=compression, significant_digits=6)
	
	# Create variables for the impurity transport results.
	imp_dens = impgrp.createVariable("imp_dens", "f4", ("gkyl_fdim", 
		"gkyl_xdim", "gkyl_ydim", "gkyl_zdim"), 
		compression=compression, significant_digits=6)
	imp_dt = impgrp.createVariable("imp_dt", "f4", 
		compression=compression)
	imp_avg_time = impgrp.createVariable("avg_time_followed", "f4", 
		compression=compression)
		
	# Create variables for the input file. Compression apparently is
	# not allowed with strings. 
	imp_mass = inputgrp.createVariable("imp_mass", "f4", 
		compression=compression)
	imp_atom_num = inputgrp.createVariable("imp_atom_num", "i4", 
		compression=compression)
	imp_init_charge = inputgrp.createVariable("imp_init_charge", "i4",
		compression=compression)
	num_imps = inputgrp.createVariable("num_imps", "i4", 
		compression=compression)
	imp_xmin = inputgrp.createVariable("imp_xmin", "f4", 
		compression=compression)
	imp_xmax = inputgrp.createVariable("imp_xmax", "f4", 
		compression=compression)	
	imp_zstart_opt = inputgrp.createVariable("imp_zstart_opt", str, 
		("str_dim",))
	imp_scaling_fact = inputgrp.createVariable("imp_scaling_fact", "f4",
		compression=compression)
	gkyl_dir = inputgrp.createVariable("gkyl_dir", str, ("str_dim",))
	gkyl_name = inputgrp.createVariable("gkyl_name", str, ("str_dim",))
	gkyl_fstart = inputgrp.createVariable("gkyl_fstart", "i4", 
		compression=compression)
	gkyl_fend = inputgrp.createVariable("gkyl_fend", "i4", 
		compression=compression)
	case_name = inputgrp.createVariable("case_name", str, ("str_dim",))
	
	# Add units to the variables.
	gkyl_x.units = "m"
	gkyl_y.units = "m"
	gkyl_z.units = "m"
	gkyl_ne.units = "m-3"
	gkyl_te.units = "eV"
	gkyl_ni.units = "m-3"
	gkyl_ti.units = "eV"
	gkyl_ex.units = "V/m"
	gkyl_ey.units = "V/m"
	gkyl_ez.units = "V/m"
	gkyl_b.units = "T"
	imp_dens.units = "m-3"
	imp_mass.units = "amu"
	imp_xmin.units = "m"
	imp_xmax.units = "m"
	imp_scaling_fact.units = "s-1"
	imp_dt.units = "s"
	imp_avg_time.units = "s"
	
	# Put variables into NetCDF files. Note strings needs to be 
	# converted to numpy array first, just the way netcdf works. 
	print("NetCDF: Saving data...")
	gkyl_x[:] = gkyl_bkg["x"]
	gkyl_y[:] = gkyl_bkg["y"]
	gkyl_z[:] = gkyl_bkg["z"]
	gkyl_ne[:] = gkyl_bkg["ne"]
	gkyl_te[:] = gkyl_bkg["te"]
	gkyl_ne[:] = gkyl_bkg["ni"]
	gkyl_te[:] = gkyl_bkg["ti"]
	gkyl_ex[:] = gkyl_bkg["elecx"]
	gkyl_ey[:] = gkyl_bkg["elecy"]
	gkyl_ez[:] = gkyl_bkg["elecz"]
	gkyl_b[:] = gkyl_bkg["b"]
	imp_dens[:] = imp_dict["imp_dens_arr"]
	imp_dt[:] = imp_dict["dt"]
	imp_avg_time[:] = imp_dict["avg_time_followed"]
	imp_mass[:] = input_opts["imp_mass"]
	imp_atom_num[:] = input_opts["imp_atom_num"]
	imp_init_charge[:] = input_opts["imp_init_charge"]
	num_imps[:] = input_opts["num_imps"]
	imp_xmin[:] = input_opts["imp_xmin"]
	imp_xmax[:] = input_opts["imp_xmax"]
	imp_zstart_opt[:] = np.array(input_opts["imp_zstart_opt"], 
		dtype="object")
	imp_scaling_fact[:] = input_opts["imp_scaling_fact"]
	gkyl_dir[:] = np.array(input_opts["gkyl_dir"])
	gkyl_name[:] = np.array(input_opts["gkyl_name"])
	gkyl_fstart[:] = input_opts["gkyl_fstart"]
	gkyl_fend[:] = input_opts["gkyl_fend"]
	case_name[:] = np.array(input_opts["case_name"])
	
	print("NetCDF: Finished.")
	
