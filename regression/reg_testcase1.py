# Regression test case #1
import pickle
import os


input_opts = {

	# Mass of impurity ion
	imp_mass        : 183.84,
	
	# Charge of impurity ion
	imp_z           : 74,
	
	# Number of impurity ions to track
	num_imps        : 1000

}


# DO NOT EDIT.
if __name__ == "__main__":
	input_fname = "{}.pickle".format(os.path.basename(__file__))
	with open(input_fname, "wb") as f:
		pickle.dump(input_opts)
