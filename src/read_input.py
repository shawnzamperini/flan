# Read input file.
import dill 
import subprocess


def read_input(input_fname):
	print("in read_input")
	
	# Run input file that shares the same name as input_fname, which 
	# will create a dilled dictionary we can store all the input
	# options in.
	input_file = "{}.py".format(input_fname)
	#exec(open(input_file).read())
	subprocess.run(["python", input_file])
	with open("{}.dill".format(input_fname), "rb") as f:
		input_opts = dill.load(f)
	return input_opts
