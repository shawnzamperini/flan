import matplotlib.pyplot as plt
import netCDF4
import numpy as np
from matplotlib.colors import Normalize, LogNorm


# Some global variables
g_fontsize = 14

class FlanPlots:

	def __init__(self, nc_path):
		
		self.nc = netCDF4.Dataset(nc_path)
	
	def load_cell_centers(self):
		"""
		Helper function to return the x, y, z cell centers. These are the 
		coordinates used to plot most data.

		Inputs: None

		Outputs:
		 - x, y, z cell centers

		Example Usage:
		 fp = flan_plots.FlanPlots("file.nc")
		 x, y, z = fp.load_cell_centers()
		"""

		x = self.nc["x"][:].data
		y = self.nc["y"][:].data
		z = self.nc["z"][:].data
		return x, y, z
	
	def load_data_frame(self, data_name, frame):
		"""
		Helper function to load a 3D set of data from the 4D data given by
		data_name at a given frame.

		Inputs:
		 - data_name (str) : Variable name of 4D data from netCDF file. 
		 - frame (int)     : The time frame to plot data for.
		
		Outputs:

		Example Usage:
		"""
		
		# Valid options that can be loaded.
		valid_opts = ["elec_dens", "elec_temp", "ion_temp", "plasma_pot", 
			"bmag"]

		# Throw an error if a valid option wasn't chosen.
		if (data_name not in valid_opts):

			# Assemble error message with the valid options.
			msg = "{} is not a valid option for ".format(data_name) + \
				"data_name. Valid options are:\n"
			for opt in valid_opts:
				msg += " {}\n".format(opt)
			raise ValueError(msg)

		return self.nc[data_name][frame].data

	def closest_index(self, arr, val):
		"""
		Helper function to return the index in arr closest to val.

		Inputs:
		 - arr (np.array) : Array of values
		 - val (float)    : Value to find nearest index for

		Outputs: Returns an int index of the closest value in arr.

		Example Usage:
		 arr = np.array([1, 2, 3, 4])
		 idx = self.closest_value(arr, 2.2)  # idx = 1
		"""

		# If arr is not a np.array make it so.
		arr = np.array(arr)

		return np.argmin(np.abs(arr - val))
	
	def get_norm(self, data, norm_type, vmin=None, vmax=None):
		"""
		Helper function to return an object to be passed into the 'norm'
		parameter in 2D matplotlib commands.

		Inputs:

		Outputs:

		Example Usage:
		"""
		
		# If vmin and/or vmax are not given, assign respectively to min/max.
		if vmin == None:
			vmin = data.min()
		if vmax == None:
			vmax = data.max()

		# Create normalization object.
		if norm_type == "linear":
			norm = Normalize(vmin=vmin, vmax=vmax)
		elif norm_type == "log":
			norm = LogNorm(vmin=vmin, vmax=vmax)

		return norm

	def plot_frame_xy(self, data_name, frame, z0, showplot=True, 
		cmap="inferno", norm_type="linear"):
		"""
		Plot data for a given frame at z=z0 in the x, y plane. data_name is
		chosen from the netCDF file, and must be 4D data (t, x, y, z). If z=z0
		is not found, the nearest value is chosen instead.

		Inputs:
		 - data_name (str) : Variable name of 4D data from netCDF file. 
		 - frame (int)     : The time frame to plot data for.
		 - z0 (float)      : The location where x, y data is plotted at.
		 - showplot (bool) : Boolean to show a plot or not.
		 - cmap (str)      : A matplotlib colormap name for the plot.
		 - norm_type (str) : The normalization of the colorbar. Options are
		                      "linear" or "log".

		Outputs: Returns a tuple of the x and y grid and the data at each
		 location.
	
		Example Usage:
		 # Returns data for user to make own plots.
		 X, Y, data_xy = self.plot_frame_xy("elec_dens", 2, 0.0, showplot=True)

		"""

		# Load cell centers.
		x, y, z = self.load_cell_centers()

		# Load data at frame.
		try:
			data = self.load_data_frame(data_name, frame)
		except ValueError as e:
			print(e)
			return None

		# Find closest z index to input z value
		z_idx = self.closest_index(z, z0) 
		print(z_idx)

		# Index data at this z location.
		data_xy = data[:,:,z_idx]

		# Grid the x, y data
		X, Y = np.meshgrid(x, y)

		# Optionally plot the data.
		if showplot:
			
			fig, ax1 = plt.subplots()

			# Plot according to the normalization option passed in.
			norm = self.get_norm(data_xy, norm_type)
			mesh = ax1.pcolormesh(X, Y, data_xy.T, cmap=cmap, norm=norm)
			
			cbar = fig.colorbar(mesh, ax=ax1)
			ax1.set_xlabel("x (m)", fontsize=g_fontsize)
			ax1.set_ylabel("y (m)", fontsize=g_fontsize)
			cbar.set_label(data_name, fontsize=g_fontsize)
			fig.tight_layout()
			fig.show()

		return X, Y, data_xy
