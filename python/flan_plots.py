import matplotlib.pyplot as plt
import netCDF4
import numpy as np
from matplotlib.colors import Normalize, LogNorm
from matplotlib import animation
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.animation import PillowWriter

# Some global variables
g_fontsize = 14

# Turn off interactive mode. This is to avoid plotting each individual plot
# when we are generating an animation.
plt.ioff()

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
		Helper function to load data_name at a given frame (if applicable, 3D
		data is returned ignoring frame).

		Inputs:
		 - data_name (str) : Variable name of data from netCDF file. 
		 - frame (int)     : The time frame to plot data for. Ignored for 3D 
		                       data.
		
		Outputs: 
		  - A 3D numpy array of x,y,z dimensions.

		Example Usage:
			# Load impurity density at frame 10
			nz = fp.load_data_frame("imp_dens", 10)
		"""

		# 3D vectors don't have a time index
		if data_name in ["b_x", "b_y", "b_z", "jacobian", "gij_00", "gij_01",
			"gij_02", "gij_11", "gij_12", "gij_22"]:
			return self.nc[data_name][:].data

		# 4D vector, index frame
		else:
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
		cmap="inferno", norm_type="linear", vmin=None, vmax=None, 
		xlabel="x (m)", ylabel="y (m)", cbar_label=None, rsep=0.0):
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

		# Time at frame, starting it at zero
		time = self.nc["time"][frame] - self.nc["time"][0]

		# Find closest z index to input z value
		z_idx = self.closest_index(z, z0) 

		# Index data at this z location.
		data_xy = data[:,:,z_idx]

		# Determine colorbar limits.
		if vmin is None: 
			vmin = data_xy.min()
		if vmax is None: 
			vmax = data_xy.max()

		# Grid the x, y data
		X, Y = np.meshgrid(x, y)
		
		# Optionally plot the data.
		if showplot:
			fig, ax1 = plt.subplots()

			# Plot according to the normalization option passed in.
			norm = self.get_norm(data_xy, norm_type, vmin=vmin, vmax=vmax)
			mesh = ax1.pcolormesh(X-rsep, Y, data_xy.T, cmap=cmap, norm=norm)
			
			cbar = fig.colorbar(mesh, ax=ax1)
			ax1.set_facecolor("grey")
			ax1.set_aspect("equal")
			ax1.set_title("{:.2f} us".format(time * 1e6))
			ax1.set_xlabel(xlabel, fontsize=g_fontsize)
			ax1.set_ylabel(ylabel, fontsize=g_fontsize)
			if cbar_label is None:
				cbar.set_label(data_name, fontsize=g_fontsize)
			else:
				cbar.set_label(cbar_label, fontsize=g_fontsize)

			fig.tight_layout()
			fig.show()

		return X, Y, data_xy

	def plot_frame_xz(self, data_name, frame, y0, showplot=True, 
		cmap="inferno", norm_type="linear", vmin=None, vmax=None, 
		xlabel="x (m)", ylabel="z (m)", cbar_label=None):
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

		# Time at frame, starting it at zero
		time = self.nc["time"][frame] - self.nc["time"][0]

		# Find closest z index to input z value
		y_idx = self.closest_index(y, y0) 

		# Index data at this z location.
		data_xz = data[:,y_idx,:]

		# Determine colorbar limits.
		if vmin is None: 
			vmin = data_xz.min()
		if vmax is None: 
			vmax = data_xz.max()

		# Grid the x, y data
		X, Z = np.meshgrid(x, z)
		
		# Optionally plot the data.
		if showplot:
			fig, ax1 = plt.subplots()

			# Plot according to the normalization option passed in.
			norm = self.get_norm(data_xz, norm_type, vmin=vmin, vmax=vmax)
			mesh = ax1.pcolormesh(X, Z, data_xz.T, cmap=cmap, norm=norm)
			
			cbar = fig.colorbar(mesh, ax=ax1)
			ax1.set_facecolor("grey")
			#ax1.set_aspect("equal")
			ax1.set_title("{:.2f} us".format(time * 1e6))
			ax1.set_xlabel(xlabel, fontsize=g_fontsize)
			ax1.set_ylabel(ylabel, fontsize=g_fontsize)
			if cbar_label is None:
				cbar.set_label(data_name, fontsize=g_fontsize)
			else:
				cbar.set_label(cbar_label, fontsize=g_fontsize)

			fig.tight_layout()
			fig.show()

		return X, Z, data_xz

	def plot_frames_xy(self, data_name, frame_start, frame_end, z0, 
		showplot=True, cmap="inferno", norm_type="linear", animate_cbar=False,
		vmin=None, vmax=None, save_path=None, xlabel="x (m)", ylabel="y (m)",
		cbar_label=None, rsep=0.0):
		"""
		Combine multiple plots from plot_frame_xy into an animation.
		"""

		# If we do not want to animate the colorbar, then we need to loop 
		# through the data ahead of time to find the limits for the colorbar.
		if not animate_cbar:
			skip_vmin = False
			skip_vmax = False
			for f in range(frame_start, frame_end+1):
				
				try:
					X, Y, data_xy = self.plot_frame_xy(data_name, f, z0, 
						showplot=False)

				# TypeError will happen if we put in an invalid option for
				# data_name.
				except TypeError:
					return None

				# Need to set vmin, vmax to values so we can use <,> on them.
				# If users pass in values for either, those take precedence and
				# we don't reassign them.
				if (f == frame_start):
					if vmin is None: 
						vmin = data_xy.min()
					else:
						skip_vmin = True
					if vmax is None: 
						vmax = data_xy.max()
					else:
						skip_vmax = True
				else:
					if (not skip_vmin and data_xy.min() < vmin): 
						vmin = data_xy.min()
					if (not skip_vmax and data_xy.max() > vmax): 
						vmax = data_xy.max()

		# Setup the first frame.
		X, Y, data_xy = self.plot_frame_xy(data_name, frame_start, z0, 
			showplot=False)

		# These commands are copied from plot_xy_frame. I wanted to have
		# plot_xy_frame return the fig and associated objects, but I couldn't
		# get around it showing the first figure in addition to the animation
		# that is later also shown. This is an undesirable effect, so we brute
		# force it by just copying the same code and reproducing it here.
		fig, ax1 = plt.subplots()
		div = make_axes_locatable(ax1)
		cax = div.append_axes("right", "5%", "5%")
		norm = self.get_norm(data_xy, norm_type, vmin=vmin, vmax=vmax)
		mesh = ax1.pcolormesh(X-rsep, Y, data_xy.T, cmap=cmap, norm=norm)
		cbar = fig.colorbar(mesh, cax=cax)
		ax1.set_facecolor("grey")
		ax1.set_aspect("equal")
		ax1.set_xlabel(xlabel, fontsize=g_fontsize)
		ax1.set_ylabel(ylabel, fontsize=g_fontsize)
		ax1.set_title("Frame {}".format(frame_start), fontsize=g_fontsize)
		if cbar_label is None:
			cbar.set_label(data_name, fontsize=g_fontsize)
		else:
			cbar.set_label(cbar_label, fontsize=g_fontsize)
		fig.tight_layout()
		fig.show()

		# Define animation function
		def animate(i):
			
			# Call for the next frame. i is being passed in zero-indexed, so to
			# get the frame we offset it from the start_frame.
			X, Y, data_xy = self.plot_frame_xy(data_name, frame_start + i, z0, 
				showplot=False)

			# Time at frame, starting it at zero
			time = self.nc["time"][i] - self.nc["time"][0]
				
			norm = self.get_norm(data_xy, norm_type, vmin=vmin, vmax=vmax)
			ax1.clear()
			mesh = ax1.pcolormesh(X-rsep, Y, data_xy.T, cmap=cmap, norm=norm)
			ax1.set_title("{:.2f} us".format(time * 1e6))
			ax1.set_xlabel(xlabel, fontsize=g_fontsize)
			ax1.set_ylabel(ylabel, fontsize=g_fontsize)
			cax.cla()
			fig.colorbar(mesh, cax=cax)
			if cbar_label is None:
				cbar.set_label(data_name, fontsize=g_fontsize)
			else:
				cbar.set_label(cbar_label, fontsize=g_fontsize)

			# Update the plot with the data from this frame.
			#ax1.clear()
			#mesh.set_array(data_xy.T)
			#vmin = data_xy.min()
			#vmax = data_xy.max()
			#mesh.set_clim(vmin, vmax)
			#cbar.update_normal(mesh)

			return mesh,

		# Generate animation
		anim = animation.FuncAnimation(fig, animate, frames=frame_end-frame_start)

		# If save_path is provided, save the file there as a gif. 
		if (save_path is not None):
			print("Saving animation...")
			anim.save(save_path + ".gif", writer=PillowWriter(fps=15))

			# Seems to be issues with this. Issue for another day I guess.
			#writer = animation.FFMpegWriter(fps=10, extra_args=["-vf", 
			#	"scale=trunc(iw/2)*2:trunc(ih/2)*2"])
			#anim.save(save_path + ".mp4", writer=writer) 

		plt.show()

	def plot_frames_xz(self, data_name, frame_start, frame_end, y0, 
		showplot=True, cmap="inferno", norm_type="linear", animate_cbar=False,
		vmin=None, vmax=None, save_path=None, xlabel="x (m)", ylabel="z (m)",
		cbar_label=None, rsep=0.0):
		"""
		Combine multiple plots from plot_frame_xz into an animation.
		"""

		# If we do not want to animate the colorbar, then we need to loop 
		# through the data ahead of time to find the limits for the colorbar.
		if not animate_cbar:
			skip_vmin = False
			skip_vmax = False
			for f in range(frame_start, frame_end+1):
				
				try:
					X, Z, data_xz = self.plot_frame_xz(data_name, f, y0, 
						showplot=False)

				# TypeError will happen if we put in an invalid option for
				# data_name.
				except TypeError:
					return None

				# Need to set vmin, vmax to values so we can use <,> on them.
				# If users pass in values for either, those take precedence and
				# we don't reassign them.
				if (f == frame_start):
					if vmin is None: 
						vmin = data_xz.min()
					else:
						skip_vmin = True
					if vmax is None: 
						vmax = data_xz.max()
					else:
						skip_vmax = True
				else:
					if (not skip_vmin and data_xz.min() < vmin): 
						vmin = data_xz.min()
					if (not skip_vmax and data_xz.max() > vmax): 
						vmax = data_xz.max()

		# Setup the first frame.
		X, Z, data_xz = self.plot_frame_xz(data_name, frame_start, y0, 
			showplot=False)

		# These commands are copied from plot_xy_frame. I wanted to have
		# plot_xy_frame return the fig and associated objects, but I couldn't
		# get around it showing the first figure in addition to the animation
		# that is later also shown. This is an undesirable effect, so we brute
		# force it by just copying the same code and reproducing it here.
		fig, ax1 = plt.subplots()
		div = make_axes_locatable(ax1)
		cax = div.append_axes("right", "5%", "5%")
		norm = self.get_norm(data_xz, norm_type, vmin=vmin, vmax=vmax)
		mesh = ax1.pcolormesh(X, Z, data_xz.T, cmap=cmap, norm=norm)
		cbar = fig.colorbar(mesh, cax=cax)
		ax1.set_facecolor("grey")
		#ax1.set_aspect("equal")
		ax1.set_xlabel(xlabel, fontsize=g_fontsize)
		ax1.set_ylabel(ylabel, fontsize=g_fontsize)
		ax1.set_title("Frame {}".format(frame_start), fontsize=g_fontsize)
		if cbar_label is None:
			cbar.set_label(data_name, fontsize=g_fontsize)
		else:
			cbar.set_label(cbar_label, fontsize=g_fontsize)
		fig.tight_layout()
		fig.show()

		# Define animation function
		def animate(i):
			
			# Call for the next frame. i is being passed in zero-indexed, so to
			# get the frame we offset it from the start_frame.
			X, Z, data_xz = self.plot_frame_xz(data_name, frame_start + i, y0, 
				showplot=False)

			# Time at frame, starting it at zero
			time = self.nc["time"][i] - self.nc["time"][0]
				
			norm = self.get_norm(data_xz, norm_type, vmin=vmin, vmax=vmax)
			ax1.clear()
			mesh = ax1.pcolormesh(X, Z, data_xz.T, cmap=cmap, norm=norm)
			ax1.set_title("{:.2f} us".format(time * 1e6))
			ax1.set_xlabel(xlabel, fontsize=g_fontsize)
			ax1.set_ylabel(ylabel, fontsize=g_fontsize)
			cax.cla()
			fig.colorbar(mesh, cax=cax)
			if cbar_label is None:
				cbar.set_label(data_name, fontsize=g_fontsize)
			else:
				cbar.set_label(cbar_label, fontsize=g_fontsize)

			# Update the plot with the data from this frame.
			#ax1.clear()
			#mesh.set_array(data_xy.T)
			#vmin = data_xy.min()
			#vmax = data_xy.max()
			#mesh.set_clim(vmin, vmax)
			#cbar.update_normal(mesh)

			return mesh,

		# Generate animation
		anim = animation.FuncAnimation(fig, animate, frames=frame_end-frame_start)

		# If save_path is provided, save the file there as a gif. 
		if (save_path is not None):
			print("Saving animation...")
			anim.save(save_path + ".gif", writer=PillowWriter(fps=15))

			# Seems to be issues with this. Issue for another day I guess.
			#writer = animation.FFMpegWriter(fps=10, extra_args=["-vf", 
			#	"scale=trunc(iw/2)*2:trunc(ih/2)*2"])
			#anim.save(save_path + ".mp4", writer=writer) 

		plt.show()

	def calc_exb_drift(self):
		"""
		Calculate and return the X,Y,Z components of the ExB drift
		"""
		pass
	
	def calc_gradb_drift(self):
		"""

		"""
		pass
	
	def calc_polarization_drift(self):
		"""

		"""
		pass
	
	def calc_curve_drift(self):
		"""

		"""
		pass
