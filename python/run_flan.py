# Top level run script for Flan. 
import sys
import read_input
import read_gkyl
import impurity_transport
import save_netcdf


def run_flan(input_fname):
	
	# Read input file. 
	input_opts = read_input.read_input(input_fname)
	
	# Read Gkeyll background data.
	gkyl_bkg = read_gkyl.read_gkyl(input_opts)
	
	# Interpolate frames if indicated.
	# To-do...
	
	# Perform impurity transport simulation. This is where the bulk of
	# the time is spent.
	imp_dict = impurity_transport.follow_impurities(input_opts, 
		gkyl_bkg)

	# Save results to a NetCDF file.
	save_netcdf.save_netcdf(input_opts, gkyl_bkg, imp_dict)

if __name__ == "__main__":
	
	if len(sys.argv) != 2:
		print("Error: Incorrect usage. Example usage:")
		print("  run_flan.py my_case")
		sys.exit()
	
	run_flan(sys.argv[1])
