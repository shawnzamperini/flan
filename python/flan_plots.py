import matplotlib.pyplot as plt
import netCDF4
import numpy as np
from matplotlib.colors import Normalize, LogNorm
from matplotlib.collections import PolyCollection
from matplotlib import animation
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.animation import PillowWriter

# Some global variables
g_fontsize = 14
amu_to_kg = 1.66e-27
elec = -1.602e-19

# Turn off interactive mode. This is to avoid plotting each individual plot
# when we are generating an animation.
plt.ioff()

class FlanPlots:

	def __init__(self, nc_path):
		
		self.nc = netCDF4.Dataset(nc_path)
		self._closed = False
	
	def close(self): 
		"""
		Close NetCDF file before leaving
		"""
		if not self._closed: 
			self.nc.close() 
			self._closed = True

	def __enter__(self): 
		return self 
	
	def __exit__(self, exc_type, exc, tb): 
		self.close()
	
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

		x = self.nc["geometry"]["x"][:].data
		y = self.nc["geometry"]["y"][:].data
		z = self.nc["geometry"]["z"][:].data
		return x, y, z
	
	def load_data_frame(self, data_name, frame, charge=1):
		"""
		Helper function to load data_name at a given frame (if applicable, 3D
		data is returned ignoring frame).

		Inputs:
		 - data_name (str) : Variable name of data from netCDF file. 
		 - frame (int)	   : The time frame to plot data for. Ignored for 3D 
							   data.
		
		Outputs: 
		  - A 3D numpy array of x,y,z dimensions.

		Example Usage:
			# Load impurity density at frame 10
			nz = fp.load_data_frame("imp_dens", 10)
		"""

		# Load data, correctly selecting the netCDF group it lives in
		# Geometry data are 3D vectors and don't have a time index
		if data_name in ["b_x", "b_y", "b_z", "J", "dXdx", "dYdx", 
			"dZdx", "dXdy", "dYdy", "dZdy", "dXdz", "dYdz", "dZdz",
			"X", "Y", "Z"]:
			return self.nc["geometry"][data_name][:].data

		# Cyclotron frequency
		elif data_name == "cyclotron_frequency":
			return self.calc_cyclo_freq()[frame]

		# Radial velocity
		elif data_name == "v_R":
			return self.calc_v_R()[frame]

		# Radial/poloidal velocities (the poloidal frame directions)
		elif data_name == "v_rad":
			return self.calc_rad_pol_comp("actual", frame=frame)[0]
		elif data_name == "v_pol":
			return self.calc_rad_pol_comp("actual", frame=frame)[1]

		# ExB drift components
		elif data_name == "ExB_X":
			return self.calc_exb_drift(frame=frame)[0]
		elif data_name == "ExB_Y":
			return self.calc_exb_drift(frame=frame)[1]
		elif data_name == "ExB_Z":
			return self.calc_exb_drift(frame=frame)[2]
		elif data_name == "ExB_rad":
			return self.calc_rad_pol_comp("ExB", frame=frame)[0]
		elif data_name == "ExB_pol":
			return self.calc_rad_pol_comp("ExB", frame=frame)[1]

		# Grad-B drift components
		elif data_name == "gradB_X":
			return self.calc_gradb_drift(charge)[0][frame]
		elif data_name == "gradB_Y":
			return self.calc_gradb_drift(charge)[1][frame]
		elif data_name == "gradB_Z":
			return self.calc_gradb_drift(charge)[2][frame]

		# Polarization drift components
		elif data_name == "pol_X":
			return self.calc_polarization_drift(charge=charge, frame=frame)[0]


		# 4D vector, index frame. Just loop through the groups until we find
		# the correct one.
		else:
			for group in ["input", "output", "background"]:
				if data_name in self.nc[group].variables.keys():
					return self.nc[group][data_name][frame].data

		# If we hit here then we couldn't find the variable
		print("Error! Cannot find variable: {}".format(data_name))
		return None

	def closest_index(self, arr, val):
		"""
		Helper function to return the index in arr closest to val.

		Inputs:
		 - arr (np.array) : Array of values
		 - val (float)	  : Value to find nearest index for

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
		xlabel="x (m)", ylabel="y (m)", cbar_label=None, rsep=0.0,
		own_data=None, charge=1, xscale=1.0, yscale=1.0, aspect="auto"):
		"""
		Plot data for a given frame at z=z0 in the x, y plane. data_name is
		chosen from the netCDF file, and must be 4D data (t, x, y, z). If z=z0
		is not found, the nearest value is chosen instead.

		Inputs:
		 - data_name (str) : Variable name of 4D data from netCDF file. 
		 - frame (int)	   : The time frame to plot data for.
		 - z0 (float)	   : The location where x, y data is plotted at.
		 - showplot (bool) : Boolean to show a plot or not.
		 - cmap (str)	   : A matplotlib colormap name for the plot.
		 - norm_type (str) : The normalization of the colorbar. Options are
							  "linear" or "log".
		 - charge (int)    : Charge of the impurity, used in calculating
							  grad-B drift.

		Outputs: Returns a tuple of the x and y grid and the data at each
		 location.
	
		Example Usage:
		 # Returns data for user to make own plots.
		 X, Y, data_xy = self.plot_frame_xy("elec_dens", 2, 0.0, showplot=True)

		"""

		# Load cell centers.
		x, y, z = self.load_cell_centers()

		# Load data at frame, either according to the variable name passed in
		# or using the data passed in via own_data.
		if own_data is None:
			try:
				data = self.load_data_frame(data_name, frame, charge)
			except ValueError as e:
				print(e)
				return None
		else:
			data = own_data[frame]

		# Time at frame, starting it at zero
		time = self.nc["geometry"]["time"][frame] \
			- self.nc["geometry"]["time"][0]

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
		X, Y = np.meshgrid(x*xscale, y*yscale)
		
		# Optionally plot the data.
		if showplot:
			fig, ax1 = plt.subplots()

			# Plot according to the normalization option passed in.
			norm = self.get_norm(data_xy, norm_type, vmin=vmin, vmax=vmax)
			mesh = ax1.pcolormesh(X-rsep, Y, data_xy.T, cmap=cmap, norm=norm)
			
			cbar = fig.colorbar(mesh, ax=ax1)
			ax1.set_facecolor("grey")
			ax1.set_aspect(aspect)
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
		xlabel="x (m)", ylabel="z (m)", cbar_label=None, own_data=None,
		aspect="auto"):
		"""
		Plot data for a given frame at z=z0 in the x, y plane. data_name is
		chosen from the netCDF file, and must be 4D data (t, x, y, z). If z=z0
		is not found, the nearest value is chosen instead.

		Inputs:
		 - data_name (str) : Variable name of 4D data from netCDF file. 
		 - frame (int)	   : The time frame to plot data for.
		 - z0 (float)	   : The location where x, y data is plotted at.
		 - showplot (bool) : Boolean to show a plot or not.
		 - cmap (str)	   : A matplotlib colormap name for the plot.
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

		# Load data at frame, either according to the variable name passed in
		# or using the data passed in via own_data.
		if own_data is None:
			try:
				data = self.load_data_frame(data_name, frame)
			except ValueError as e:
				print(e)
				return None
		else:
			data = own_data[frame]

		# Time at frame, starting it at zero
		time = self.nc["geometry"]["time"][frame] \
			- self.nc["geometry"]["time"][0]

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
			ax1.set_aspect(aspect)
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
		cbar_label=None, rsep=0.0, aspect="auto", own_data=None):
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
						showplot=False, own_data=own_data)

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
			showplot=False, own_data=own_data)

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
		ax1.set_aspect(aspect)
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
				showplot=False, own_data=own_data)

			# Time at frame, starting it at zero
			time = self.nc["geometry"]["time"][i] \
				- self.nc["geometry"]["time"][0]
				
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
		cbar_label=None, rsep=0.0, own_data=None, aspect="auto"):
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
						showplot=False, own_data=own_data)

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
			showplot=False, own_data=own_data)

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
		ax1.set_aspect(aspect)
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
				showplot=False, own_data=own_data)

			# Time at frame, starting it at zero
			time = self.nc["geometry"]["time"][i] \
				- self.nc["geometry"]["time"][0]
				
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

	def calc_v_R(self):
		"""
		Calculate and return the radial velocity

		Outputs: Returns an array of the radial velocity [t,x,y,z]
		"""
		
		# Pull out some arrays
		vX = self.nc["output"]["v_X"][:]
		vY = self.nc["output"]["v_Y"][:]
		X = self.nc["geometry"]["X"][:]
		Y = self.nc["geometry"]["Y"][:]
		vR_mag = np.sqrt(np.square(vX) + np.square(vY))

		# If scalar_proj < 0 we need to multiply vR by -1
		scalar_proj = (vX * X + vY * Y)
		vR_sign = np.ones(scalar_proj.shape)
		vR_sign[scalar_proj < 0] = -1

		return vR_mag * vR_sign

	def calc_exb_drift(self, frame=None):
		"""
		Calculate and return the X,Y,Z components of the ExB drift

		Outputs: Returns a tuple of the X,Y,Z ExB drift components
		"""
		
		# Pull out some arrays, getting a specific frame if desired
		if frame == None:
			EX = self.nc["background"]["E_X"][:]
			EY = self.nc["background"]["E_Y"][:]
			EZ = self.nc["background"]["E_Z"][:]
			BX = self.nc["background"]["B_X"][:]
			BY = self.nc["background"]["B_Y"][:]
			BZ = self.nc["background"]["B_Z"][:]
		else:
			EX = self.nc["background"]["E_X"][frame]
			EY = self.nc["background"]["E_Y"][frame]
			EZ = self.nc["background"]["E_Z"][frame]
			BX = self.nc["background"]["B_X"][frame]
			BY = self.nc["background"]["B_Y"][frame]
			BZ = self.nc["background"]["B_Z"][frame]

		Bsq = np.square(BX) + np.square(BY) + np.square(BZ)

		# Perpendicular component of the electric field
		scalar_proj = (EX * BX + EY * BY + EZ * BZ)
		EparX = scalar_proj * BX / Bsq
		EparY = scalar_proj * BY / Bsq
		EparZ = scalar_proj * BZ / Bsq
		EperpX = EX - EparX
		EperpY = EY - EparY
		EperpZ = EZ - EparZ

		# Calculate each component and return 
		v_exb_X = (EperpY * BZ - EperpZ * BY) / Bsq
		v_exb_Y = -(EperpX * BZ - EperpZ * BX) / Bsq
		v_exb_Z = (EperpX * BY - EperpY * BX) / Bsq
		return v_exb_X, v_exb_Y, v_exb_Z
	
	def calc_gradb_drift(self, charge, frame=None):
		"""
		Calculate and return the X,Y,Z components of the grad-B drift

		Inputs
			- charge: Charge of the impurity

		Outputs: Returns a tuple of the X,Y,Z grad-B drift components
		"""
		
		# Pull out some arrays
		if frame == None:
			vX = self.nc["output"]["v_X"][:]
			vY = self.nc["output"]["v_Y"][:]
			vZ = self.nc["output"]["v_Z"][:]
			BX = self.nc["background"]["B_X"][:]
			BY = self.nc["background"]["B_Y"][:]
			BZ = self.nc["background"]["B_Z"][:]
			gradBX = self.nc["background"]["dBdX"][:]
			gradBY = self.nc["background"]["dBdY"][:]
			gradBZ = self.nc["background"]["dBdZ"][:]
		else:
			vX = self.nc["output"]["v_X"][frame]
			vY = self.nc["output"]["v_Y"][frame]
			vZ = self.nc["output"]["v_Z"][frame]
			BX = self.nc["background"]["B_X"][frame]
			BY = self.nc["background"]["B_Y"][frame]
			BZ = self.nc["background"]["B_Z"][frame]
			gradBX = self.nc["background"]["dBdX"][frame]
			gradBY = self.nc["background"]["dBdY"][frame]
			gradBZ = self.nc["background"]["dBdZ"][frame]

		mz_kg = self.nc["input"]["imp_mass_amu"][:] * amu_to_kg

		# Magnetic field magnitude squared
		Bsq = np.square(BX) + np.square(BY) + np.square(BZ)

		# Perpendicular to B velocity. Can obtain the perpendicular components
		# by subtracting the parallel projection of the impurity velocity from 
		# the total velocity vector.
		scalar_proj = (vX * BX + vY * BY + vZ * BZ)
		vparX = scalar_proj * BX / Bsq
		vparY = scalar_proj * BY / Bsq
		vparZ = scalar_proj * BZ / Bsq
		vperpX = vX - vparX
		vperpY = vY - vparY
		vperpZ = vZ - vparZ
		vperp_sq = np.square(vperpX) + np.square(vperpY) + np.square(vperpZ)

		# Cyclotron frequency
		omega_c = np.abs(charge * elec) * np.sqrt(Bsq) / mz_kg

		# Grad-B drift components for a *positively* charged particle. There
		# is a -/+ sign in front, where the + applies to positive charge. 
		# See Freidberg Eq. 8.105.
		v_gradB_X = vperp_sq / (2 * omega_c * Bsq) * (BY * gradBZ - BZ * gradBY)
		v_gradB_Y = -vperp_sq / (2 * omega_c * Bsq) * (BX * gradBZ - BZ * gradBX)
		v_gradB_Z = vperp_sq / (2 * omega_c * Bsq) * (BX * gradBY - BY * gradBX)
		return v_gradB_X, v_gradB_Y, v_gradB_Z
	
	def calc_polarization_drift(self, charge, frame=None):
		"""
		Calculate and return the X,Y,Z components of the polarization drift

		Inputs
			- charge: Charge of the impurity

		Outputs: Returns a tuple of the X,Y,Z polarization drift components
		"""

		# Pull out some arrays
		if frame == None:
			vX = self.nc["output"]["v_X"][:]
			vY = self.nc["output"]["v_Y"][:]
			vZ = self.nc["output"]["v_Z"][:]
			BX = self.nc["background"]["B_X"][:]
			BY = self.nc["background"]["B_Y"][:]
			BZ = self.nc["background"]["B_Z"][:]
			t = self.nc["geometry"]["time"][:]

		# If we just specify a frame, we will still need the previous frame
		# so we can do a forward difference for the acceleration calculation
		else:
			vX = self.nc["output"]["v_X"][frame-1:frame+1]
			vY = self.nc["output"]["v_Y"][frame-1:frame+1]
			vZ = self.nc["output"]["v_Z"][frame-1:frame+1]
			BX = self.nc["background"]["B_X"][frame-1:frame+1]
			BY = self.nc["background"]["B_Y"][frame-1:frame+1]
			BZ = self.nc["background"]["B_Z"][frame-1:frame+1]
			t = self.nc["geometry"]["time"][frame-1:frame+1]

		mz_kg = self.nc["input"]["imp_mass_amu"][:] * amu_to_kg

		# Magnetic field magnitude squared
		Bsq = np.square(BX) + np.square(BY) + np.square(BZ)

		# Magnetic field unit vector
		bX = BX / Bsq
		bY = BY / Bsq
		bZ = BZ / Bsq
		
		# Cyclotron frequency
		omega_c = np.abs(charge * elec) * np.sqrt(Bsq) / mz_kg

		# The polarization drift includes the acceleration, which is not 
		# actually stored but we can estimate it easily with the change in
		# average velocity over each time frame.
		aX = np.zeros(vX.shape)
		aY = np.zeros(vY.shape)
		aZ = np.zeros(vZ.shape)
		for i in range(1, len(t)):
			dt = t[i] - t[i-1]
			aX[i] = (vX[i] - vX[i-1]) / dt
			aY[i] = (vY[i] - vY[i-1]) / dt
			aZ[i] = (vZ[i] - vZ[i-1]) / dt

		# Polarization drift components for a *positively* charged particle.
		# There is a -/+ sign in front, where the + applies to positive charge. 
		# See Freidberg Eq. 8.105.
		v_pol_X = (bY * aZ - bZ * aY) / omega_c
		v_pol_Y = -(bX * aZ - bZ * aX) / omega_c
		v_pol_Z = (bX * aY - bY * aX) / omega_c
	
		# Return the components at the desired frame, if necessary, which
		# is just index = 1 since 0 is the previous frame used for calculating
		# derivatives.
		if frame == None:
			return v_pol_X, v_pol_Y, v_pol_Z
		else:
			return v_pol_X[1], v_pol_Y[1], v_pol_Z[1]
	
	def calc_curvature_drift(self, charge):
		"""
		Calculate and return the X,Y,Z components of the curvature drift

		Inputs
			- charge: Charge of the impurity

		Outputs: Returns a tuple of the X,Y,Z curvature drift components
		"""
		
		# Pull out some arrays
		vX = self.nc["output"]["v_X"][:]
		vY = self.nc["output"]["v_Y"][:]
		vZ = self.nc["output"]["v_Z"][:]
		BX = self.nc["background"]["B_X"][:]
		BY = self.nc["background"]["B_Y"][:]
		BZ = self.nc["background"]["B_Z"][:]
		mz_kg = self.nc["input"]["imp_mass_amu"][:] * amu_to_kg

		# Magnetic field magnitude squared
		Bsq = np.square(BX) + np.square(BY) + np.square(BZ)
		B = np.sqrt(Bsq)

		# Cyclotron frequency
		omega_c = np.abs(charge * elec) * np.sqrt(Bsq) / mz_kg

		# Parallel to B velocity components
		vparX = (vX * BX + vY * BY + vZ * BZ) * BX / Bsq
		vparY = (vX * BX + vY * BY + vZ * BZ) * BY / Bsq
		vparZ = (vX * BX + vY * BY + vZ * BZ) * BZ / Bsq
		vpar_sq = np.square(vparX) + np.square(vparY) + np.square(vparZ)

		# Curvature vector. This is normally defined in terms of the magnetic
		# field, but we have the X, Y coordinates so can just use that here.
		X = self.nc["geometry"]["X"][:]
		Y = self.nc["geometry"]["Y"][:]
		R_sq = np.square(X) + np.square(Y)

		# Calculate each component of curvature drift and return
		v_curv_X = np.zeros(vX.shape)
		v_curv_Y = np.zeros(vY.shape)
		v_curv_Z = np.zeros(vZ.shape)
		for i in range(vX.shape[0]):
			v_curv_X[i] = (Y * BZ[i]) * vpar_sq[i] / (omega_c[i] * R_sq * B[i])
			v_curv_Y[i] = -(X * BZ[i]) * vpar_sq[i] / (omega_c[i] 
				* R_sq[i] * B[i])
			v_curv_Z[i] = (X * BY[i] - Y * BX[i]) * vpar_sq[i] / (omega_c[i] 
				* R_sq * B[i])
		return v_curv_X, v_curv_Y, v_curv_Z
	
	def calc_cyclo_freq(self):
		"""
		Calculate average cyclotron frequency using the average charges. 
		"""

		# Pull out some arrays
		BX = self.nc["background"]["B_X"][:]
		BY = self.nc["background"]["B_Y"][:]
		BZ = self.nc["background"]["B_Z"][:]
		qz = self.nc["output"]["qz"][:]
		mz_kg = self.nc["input"]["imp_mass_amu"][:] * amu_to_kg

		# Magnetic field magnitude squared
		Bsq = np.square(BX) + np.square(BY) + np.square(BZ)

		# Cyclotron frequency
		omega_c = np.abs(qz * elec) * np.sqrt(Bsq) / mz_kg
		return omega_c

	def calc_gyrorad(self):
		"""
		Calculate average gyroradius using the average velocity and charges. 
		"""

		# Pull out some arrays
		vX = self.nc["output"]["v_X"][:]
		vY = self.nc["output"]["v_Y"][:]
		vZ = self.nc["output"]["v_Z"][:]
		BX = self.nc["background"]["B_X"][:]
		BY = self.nc["background"]["B_Y"][:]
		BZ = self.nc["background"]["B_Z"][:]
		qz = self.nc["output"]["qz"][:]
		mz_kg = self.nc["input"]["imp_mass_amu"][:] * amu_to_kg

		# Magnetic field magnitude squared
		Bsq = np.square(BX) + np.square(BY) + np.square(BZ)

		# Perpendicular to B velocity. Can obtain the perpendicular components
		# by subtracting the parallel projection of the impurity velocity from 
		# the total velocity vector.
		scalar_proj = (vX * BX + vY * BY + vZ * BZ)
		vparX = scalar_proj * BX / Bsq
		vparY = scalar_proj * BY / Bsq
		vparZ = scalar_proj * BZ / Bsq
		vperpX = vX - vparX
		vperpY = vY - vparY
		vperpZ = vZ - vparZ
		vperp = np.sqrt(np.square(vperpX) + np.square(vperpY) 
			+ np.square(vperpZ))

		# Cyclotron frequency
		omega_c = np.abs(qz * elec) * np.sqrt(Bsq) / mz_kg

		# Larmor radius
		return vperp / omega_c
	
	def calc_E_R(self):
		"""
		Calculates the radial electric field, ensuring to assign the correct 
		sign. Positive = radially outwards, nagative = radially inwards.
		"""

		# Pull out arrays
		X = self.nc["geometry"]["X"][:].data
		Y = self.nc["geometry"]["Y"][:].data
		EX = self.nc["background"]["E_X"][:].data
		EY = self.nc["background"]["E_Y"][:].data

		# Loop through each time
		ER = np.zeros(EX.shape)
		for i in range(EX.shape[0]):

			# Scalar projection of E onto R
			scalar_proj = X * EX[i] + Y * EY[i]

			# Direction scalar to determine if ER points radially outwards (+) or 
			# inwards (-)
			dir_scalar = np.ones(X.shape)
			dir_scalar[scalar_proj < 0] = -1

			# Put it all together
			ER[i] = dir_scalar * np.sqrt(np.square(EX[i]) + np.square(EY[i]))

		return ER

	# Refactor this to accept any general vector, make it lower level.
	# As written this is a little messy.
	def calc_rad_pol_comp(self, which, frame=None):
		"""
		Calculate the radial and poloidal components of a vector field. In
		this context, radial is the cross-field direction and not the
		cylindrical radial coordinate.

		which:	Which velocity to calculate radial/poloidal components
				for. One of: 
					ExB, 
					gradB, 
					polarization, 
					curvature
					actual. 

		Returns both the radial and poloidal components as either a 4D (all
		times, frame = None) or 3D (specific time, frame = integer) array. Slots
		into the other functions easily.
		"""

		# Validate input
		valid_which = ["actual", "exb", "gradb", "polarization", "curvature"]
		if which.lower() not in valid_which:
			print(f"Error! which = {which} not valid. Valid options: "
				  + ", ".join(valid_which))
			return None

		# Geometry arrays (nx, ny, nz)
		X = self.nc["geometry"]["X"][:].data
		Y = self.nc["geometry"]["Y"][:].data
		Z = self.nc["geometry"]["Z"][:].data
		R = np.sqrt(X**2 + Y**2)

		# Magnetic field arrays
		if frame is None:
			BX = self.nc["background"]["B_X"][:]   # (nt, nx, ny, nz)
			BY = self.nc["background"]["B_Y"][:]
			BZ = self.nc["background"]["B_Z"][:]
		else:
			BX = np.expand_dims(self.nc["background"]["B_X"][frame], axis=0)
			BY = np.expand_dims(self.nc["background"]["B_Y"][frame], axis=0)
			BZ = np.expand_dims(self.nc["background"]["B_Z"][frame], axis=0)

		# Velocity components
		if which.lower() == "actual":
			if frame is None:
				vX = self.nc["output"]["v_X"][:]   # (nt, nx, ny, nz)
				vY = self.nc["output"]["v_Y"][:]
				vZ = self.nc["output"]["v_Z"][:]
			else:
				vX = np.expand_dims(self.nc["output"]["v_X"][frame], axis=0)
				vY = np.expand_dims(self.nc["output"]["v_Y"][frame], axis=0)
				vZ = np.expand_dims(self.nc["output"]["v_Z"][frame], axis=0)

		elif which.lower() == "exb":
			vX, vY, vZ = calc_exb_drift(self, frame=frame)

		# Stack into vector fields
		B = np.stack([BX, BY, BZ], axis=-1)		 # (nt, nx, ny, nz, 3)
		v = np.stack([vX, vY, vZ], axis=-1)

		# Toroidal unit vector e_phi = (-y, x, 0) / R
		e_phi = np.zeros(X.shape + (3,))
		e_phi[..., 0] = -Y / R
		e_phi[..., 1] =  X / R
		e_phi[..., 2] =  0.0

		# Normalize
		e_phi /= np.linalg.norm(e_phi, axis=-1, keepdims=True)

		# Broadcast e_phi to match time dimension. Broadcast requires the
		# rightmost index match, which in this case is 3 for both arrays. Now
		# e_phi behaves like an (nt, nx, ny, nz, 3) shaped array without 
		# actually using any additional memory.
		e_phi = np.broadcast_to(e_phi, B.shape)

		# Poloidal direction = B × e_phi
		e_pol = np.cross(B, e_phi)
		e_pol /= np.linalg.norm(e_pol, axis=-1, keepdims=True)
		
		# Radial direction = e_phi × e_pol
		e_r = np.cross(e_phi, e_pol)
		e_r /= np.linalg.norm(e_r, axis=-1, keepdims=True)

		# Project velocity onto radial and poloidal directions. This is just
		# the dot product. v*e_r is the 3 terms in the dot product, then we
		# sum them along the coordinate axes to finish the dot product.
		v_rad = np.sum(v * e_r, axis=-1)
		v_pol = np.sum(v * e_pol, axis=-1)

		# If a single frame was requested, return 3D arrays
		if frame is not None:
			return v_rad[0], v_pol[0]

		# Otherwise return full 4D arrays
		return v_rad, v_pol		
	
	def validate_drift_test(self):
		"""
		"""
		
		from scipy.signal import find_peaks

		# Check that particle tracks were on
		if self.nc["input"]["save_track"][:] != "on":
			print("Error! Particle track was not saved. Rerun \
			with save_track =\" on\"")
			return None

		# Load x,y,z versus t
		t = self.nc["output"]["track_t"][:].data
		x = self.nc["output"]["track_x"][:].data
		y = self.nc["output"]["track_y"][:].data
		z = self.nc["output"]["track_z"][:].data

		# Load particle traits
		q = self.nc["input"]["imp_init_charge"][0] * (-elec)
		m_amu = self.nc["input"]["imp_mass_amu"][0]
		m = m_amu * amu_to_kg

		# Constants for plotting
		fontsize = 16
		figsize = (5, 4)

		# Check which test was performed. For reference:
		# 0 = Simple gyration test to show a particle gyrates in a constant
		#     magnetic field
		# 1 = Test to show a gyrating particle will ExB drift in a constant
		#     electric and magnetic field
		# 2 = Test to show a gyrating particle will grad-B drift in a magentic
		#     field with a linear gradient
		# 3 = Test to show a gyrating particle will have a polarization drift
		#     in a time-varying electric field and constant magnetic field
		# 4 = Test to show a gyrating particle will have a curvature drift
		#     in a toroidal magnetic field
		# 5 = Test to show that the fluid friction force manifests  from the
		#     collision model
		if self.nc["input"]["test_opt"][0] == "gyrate":
			pass

		# Input file: tests/test_exb.cpp
		# Background plasma is constant BZ = 1 and EY = 500. Drift is in the
		# x direction.
		elif self.nc["input"]["test_opt"][0] == "exb":
			
			print("ExB drift test case detected")

			# Hardcoded values that this test uses
			vY0 = 25000
			EY = 500
			BZ = 1.0

			# Cyclotron period = qB / m
			omega = q * BZ / m

			# Calculate the analytic drift value
			analytic_drift = EY / BZ

			# Drift in x direction so we want to find those peaks
			drift_coord = x
			drift_ylabel = "x (m)"

		# Input file: tests/test_gradb.cpp
		elif self.nc["input"]["test_opt"][0] == "gradb":
			
			print("Grad-B drift test case detected")

			# Hardcoded values that this test uses (see read_test.cpp)
			vY0 = 25000
			EY = 0
			BZ = 1.0
			dBZdy = 5.0  # 5 T/m
			omega = q * BZ / m

			# Analytic estimate of the Grad-B drift velocity
			analytic_drift = -m * vY0**2 * dBZdy / (2 * q * BZ**2)	

			# Drift in x direction so we want to find those peaks
			drift_coord = x
			drift_ylabel = "x (m)"

		# Input file: tests/test_polarization.cpp
		elif self.nc["input"]["test_opt"][0] == "polarization":

			print("Polarization drift test case detected")

			# Hardcoded values that this test uses (see read_test.cpp)
			vY0 = 25000
			EY = 5000
			BZ = 1.0
			dEYdt = -10000000  # 1 V/us
			omega = q * BZ / m

			# Analytic estimate of the polarization drift velocity
			analytic_drift = m * dEYdt / (q * BZ**2)

			# Drift in y direction so we want to find those peaks
			drift_coord = y
			drift_ylabel = "y (m)"

		# Input file: tests/test_curvature.cpp
		elif self.nc["input"]["test_opt"][0] == "curvature":

			print("Curvature drift test case detected")

			# Hardcoded values that this test uses (see read_test.cpp)
			B0 = 10.0
			v_par = 5000
			R = x.mean()
			dt = self.nc["input"]["imp_time_step"][:].data[0]

			omega = q * B0 / m

			# Drift in y=Z direction so we want to find those peaks
			drift_coord = y
			drift_ylabel = "Z (m)"

			# Analytic estimate of the curvature drift velocity
			analytic_drift = v_par**2 / omega / R

			# This test gets an additional plot to demonstrate the numerical
			# radial drift that can appear in curved geometries
			numerical_drift = -v_par**2 * dt / R

			peaks, properties = find_peaks(x)
			t_peaks = t[peaks]
			drift_peaks = x[peaks]
			v_drift, b = np.polyfit(t_peaks, drift_peaks, 1)
			tfit = np.linspace(t.min(), t.max(), 100)
			drift_fit = v_drift * tfit + b

			fig2, ax2 = plt.subplots(figsize=figsize)
			ax2.set_title("Numerical Drift", fontsize=fontsize)
			ax2.set_xlabel("Time (s)", fontsize=fontsize)
			ax2.set_ylabel("R (m)", fontsize=fontsize)
			ax2.plot(t, x, color="k", lw=3, label="Flan")
			ax2.plot(tfit, drift_fit, color="k", lw=3, linestyle="--", 
				label="Analytic")
			ax2.legend()
			fig2.tight_layout()
			fig2.show()

			abs_err = abs(v_drift - numerical_drift)
			rel_err = abs_err / abs(numerical_drift)
			print(f"Cartesian Numerical Drift\n"
				f"  Computed : {v_drift: .2f} m/s\n"
				f"  Analytic : {analytic_drift: .2f} m/s\n"
				f"  Abs. err : {abs_err: .3f} m/s\n"
				f"  Rel. err : {rel_err: .2e}\n")

		# Slope from linear fit to peaks is simulated drift
		peaks, properties = find_peaks(drift_coord)
		t_peaks = t[peaks]
		drift_peaks = drift_coord[peaks]
		v_drift, b = np.polyfit(t_peaks, drift_peaks, 1)
		tfit = np.linspace(t.min(), t.max(), 100)
		drift_fit = v_drift * tfit + b

		fig, ax1 = plt.subplots(figsize=figsize)
		ax1.set_xlabel("Time (s)", fontsize=fontsize)
		ax1.set_ylabel(drift_ylabel, fontsize=fontsize)
		ax1.plot(t, drift_coord, color="k", lw=3, label="Flan")
		ax1.plot(tfit, drift_fit, color="k", lw=3, linestyle="--", 
			label="Analytic")
		ax1.legend()
		fig.tight_layout()
		fig.show()

		# Print summary
		abs_err = abs(v_drift - analytic_drift)
		rel_err = abs_err / abs(analytic_drift)

		print(f"\n"
			f"  Computed : {v_drift: .2f} m/s\n"
			f"  Analytic : {analytic_drift: .2f} m/s\n"
			f"  Abs. err : {abs_err: .3f} m/s\n"
			f"  Rel. err : {rel_err: .2e}\n")
	
	def validate_coll_test(self):
		"""
		"""

	
	def validate_test(self, show_nanbu_s=True):
		"""
		Make a plots of the hardcoded tests cases within Flan that ensures the 
		physics models within Flan are working correctly.
		"""

		from scipy.signal import find_peaks

		# Constants for plotting
		fontsize = 16
		figsize = (5, 4)

		# Some simulation constants
		q = self.nc["input"]["imp_init_charge"][0] * (-elec)
		m_amu = self.nc["input"]["imp_mass_amu"][0]
		m = m_amu * amu_to_kg

		# Check that a test was indeed performed via bkg_source
		if self.nc["input"]["bkg_source"][:] != "test":
			print("Error! This is not a test simulation.")
			return None

		# Call drift test routine
		if self.nc["input"]["test_opt"][0] in ["gyrate", "exb", "gradb",
			"polarization", "curvature"]:
			self.validate_drift_test()

		# Input file: tests/test_friction_force.cpp
		elif self.nc["input"]["test_opt"][0] == "friction_force":

			# Hardcoded values that this test uses (see read_test.cpp)
			BX = 1.0
			uX = 1000
			Ti = 1
			mi = 2.014
			ni = 1e20
			#ln_alpha = 10
			ln_alpha = 3.5 # n=1e20, T=1, u=1000
			Z = self.nc["input"]["imp_init_charge"][0] # Ion/recomb is OFF

			# Some additional values and arrays
			t = self.nc["geometry"]["time"][:]
			#t0 = t[0]
			#vx = self.nc["output"]["track_vx"][:]

			# Average velocity of particles that all start at the same
			# time and place. Ignore non-zero locations.
			vx = np.zeros(len(t))
			vx_std = np.zeros(len(t))
			s = np.zeros(len(t))
			for i in range(len(t)):
				with_counts = self.nc["output"]["Nz"][i] > 0
				vx[i] = self.nc["output"]["v_X"][i][with_counts].mean()
				vx_std[i] = self.nc["output"]["v_X"][i][with_counts].std()
				if show_nanbu_s:
					s[i] = self.nc["output"]["nanbu_s"][i][with_counts].mean()

			# Remove NaN data
			notnan = ~np.isnan(vx)
			t = t[notnan]
			vx = vx[notnan]
			vx_std = vx_std[notnan]
			s = s[notnan]

			from scipy.optimize import curve_fit

			def vx_model(t, nu):
				return uX * (1 - np.exp(-nu * t))

			# --- 3. Fit ν_eff ---
			popt, pcov = curve_fit(vx_model, t, vx, p0=[1e4])
			nu_eff = popt[0]

			print("Fitted ν_eff =", nu_eff)

			# --- 4. Plot to verify ---
			t_fit = np.linspace(t.min(), t.max(), 500)
			vx_fit = vx_model(t_fit, nu_eff)

			# Running average window filter
			dt = t[1] - t[0]
			window = int(1e-6 / dt)
			#window = 50
			#vx_avg = np.convolve(vx, np.ones(window)/window, mode='same')
			
			# Stopping power time from Stangeby 6.35, in s
			tau_s = 1.47e13 * m_amu * Ti * np.sqrt(Ti / mi) \
				/ ((1 + mi / m_amu) * ni * Z**2 * ln_alpha)
			print(f"tau_s = {tau_s:.2e} s")
			print(f"nu_s = {1/tau_s:.2e} s")
			print(f"nanbu_s = {s.mean():.2f}")

			# Friction force
			#ff_avg = m * (uX - vz_t_avg) / tau_s
			#ff_max = m * (uX - vz_t_max) / tau_s

			# Simple first-order ODE to solve for vz from F=ma
			vx_analytic = uX * (1 - np.exp(-t / tau_s))

			for i in range(len(t)):
				t_i = t[i]
				vx_i = vx[i]
				vx_std_i = vx_std[i]
				vxa_i = vx_analytic[i]
				print(f"{t_i:.2e} {vx_i:.2e} {vx_std_i:.2e} {vxa_i:.2e}") 

			fig, ax1 = plt.subplots(figsize=figsize)
			ax1.fill_between(t, vx-vx_std, vx+vx_std, color="tab:red", alpha=0.3)
			ax1.plot(t, vx, color="tab:red", lw=3, label="Flan")
			ax1.plot(t_fit, vx_fit, color="k", lw=3)
			ax1.plot(t_fit, vx_fit, color="tab:red", lw=2)
			#ax1.plot(t, vx_avg, color="tab:red", lw=3, label="Flan")
			ax1.plot(t, vx_analytic, color="k", lw=3, linestyle="--", 
				label="Analytic")
			ax1.set_xlabel("Time (s)", fontsize=fontsize)
			ax1.set_ylabel("vX (m/s)", fontsize=fontsize, color="tab:red")
			ax1.set_ylim([0, 2000])

			if show_nanbu_s:
				ax11 = ax1.twinx()
				ax11.plot(t, s, lw=3, color="k")
				ax11.plot(t, s, lw=2, color="tab:purple")
				ax11.set_ylabel("Nanbu - s", fontsize=fontsize, 
					color="tab:purple")
				ax11.set_ylim([0, 40])

			fig.tight_layout()
			fig.show()

		else:
			
			print("Error: test_opt = {} not recognized" \
				.format(self.nc["input"]["test_opt"][0]))
			return None
	
	def plot_frame_RZ(self, data_name, frame, nodes_path, charge=1, 
		gfile_path=None, cmap="inferno", norm_type="linear", 
		show_pol_ang=True, vmin=None, vmax=None):
		"""
		Plot toroidally averaged data in R, Z frame.
		"""
		pass

		# Needed to read nodes file. Built into flan conda environment.
		import postgkyl as pg

		# Pull out some arrays for easy access
		x = self.nc["geometry"]["x"][:]  # psi
		y = self.nc["geometry"]["y"][:]  # alpha
		z = self.nc["geometry"]["z"][:]  # poloidal theta, I think!

		# Average over y because we treat y as the toroidal coordinate, so this
		# is a toroidal average.
		data_yavg = self.load_data_frame(data_name, frame, charge).mean(axis=1)

		# Load gfile. This is just used to plot the flux surfaces and is not 
		# required (but obviously very useful)
		include_gfile = False
		if gfile_path is not None:
			
			from read_gfile import read_gfile

			gfile = read_gfile(gfile_path)
			rgrid = gfile["rgrid"]
			zgrid = gfile["zgrid"]
			psigrid = gfile["psi_grid"]
			psirz = gfile["psirz"].T  # Indexing Fortran --> C (which is assumed in numpy)
			R_axis = gfile["rmaxis"]
			Z_axis = gfile["zmaxis"]
			psi_axis = gfile["simag"]
			psi_lcfs = gfile["sibry"]
			qpsi = gfile["qpsi"]
			rlim = gfile["rlim"]
			zlim = gfile["zlim"]

			# Set flag so plot options execute
			include_gfile = True

		# Load nodes file.
		nodes_gdata = pg.GData(nodes_path)
		nodes = nodes_gdata.get_values()  # shape (x+1, y+1, z+1, 3)

		# nodes contains the R,Z,phi coordinates of the nodes, which are the 0,1,2
		# index of the last dimension, respectively. We average over the y 
		# (toroidal) index for the average coordinates in the R, Z plane.
		R_nodes = nodes[:, :, :, 0].mean(axis=1)
		Z_nodes = nodes[:, :, :, 1].mean(axis=1)

		# x and z node coordinates and cell centers
		x_nodes = nodes_gdata.get_grid()[0]
		z_nodes = nodes_gdata.get_grid()[2]
		x_centers = 0.5 * (x_nodes[:-1] + x_nodes[1:])
		z_centers = 0.5 * (z_nodes[:-1] + z_nodes[1:])

		Nx = len(x_nodes) - 1  # number of cells in x
		Nz = len(z_nodes) - 1  # number of cells in z

		# R_nodes and Z_nodes are already shape (Nx, Nz) = cell centers from y-averaging
		# But we need node positions at the node grid, not cell centers.
		# Use the raw nodes array directly (before averaging over y):
		# nodes shape is (Nx, Ny, Nz, 3), take y=0 (all y are identical as we confirmed)
		R_node_grid = nodes[:, 0, :, 0]  # shape (Nx, Nz) - but these are NODE positions
		Z_node_grid = nodes[:, 0, :, 1]  # shape (Nx, Nz)

		# So nodes gives R,Z at the LOW corner of each cell, not at all 4 
		# corners. We need to construct the 4 corners from adjacent cells.

		# Build quads: for each cell (i,j), the 4 corners in (R,Z) are:
		# (i,j), (i+1,j), (i+1,j+1), (i,j+1)
		# We need R,Z at all node indices including the last row/col.
		# Extend R_node_grid to (Nx+1, Nz+1) by extrapolating the boundary:

		def extend_node_grid(G):
			"""
			Extend a (Nx, Nz) node grid to (Nx+1, Nz+1) by linear extrapolation.
			"""

			# Add one row at the end (x direction)
			last_row = 2 * G[-1, :] - G[-2, :]
			G = np.vstack([G, last_row[np.newaxis, :]])

			# Add one column at the end (z direction)
			last_col = 2 * G[:, -1] - G[:, -2]
			G = np.hstack([G, last_col[:, np.newaxis]])
			return G

		R_ext = extend_node_grid(R_node_grid)  # shape (Nx+1, Nz+1)
		Z_ext = extend_node_grid(Z_node_grid)  # shape (Nx+1, Nz+1)

		# Now build the polygon list
		verts = []
		colors = []
		for i in range(Nx):
			for j in range(Nz-1):  # I think -1 is needed here

				# 4 corners, ordered as a closed polygon
				quad = np.array([
					[R_ext[i,   j  ], Z_ext[i,   j  ]],
					[R_ext[i+1, j  ], Z_ext[i+1, j  ]],
					[R_ext[i+1, j+1], Z_ext[i+1, j+1]],
					[R_ext[i,   j+1], Z_ext[i,   j+1]],
				])
				verts.append(quad)

				# The nodes are not necessarily on the same grid that the
				# simulation is carried out on, so we need to grab the nearest 
				# value. This is because a higher resolution grid allows us to
				# better resolves the shape of the flux surface in plots. So
				# we need to find nearest simulation indices to each center
				# from the nodes file.
				i_sim = np.argmin(np.abs(x - x_centers[i]))
				j_sim = np.argmin(np.abs(z - z_centers[j]))
				colors.append(data_yavg[i_sim, j_sim])

		colors = np.array(colors)

		# Mask out non-positive values for log scale
		if norm_type == "log":
			valid = colors > 0
			verts_valid = [v for v, ok in zip(verts, valid) if ok]
			colors_valid = colors[valid]
		else:
			verts_valid = verts
			colors_valid = colors

		# Assign vmix/vamx
		if vmin is None and norm_type == "log":
			vmin = np.nanmin(colors_valid)	
		if vmax is None and norm_type == "log":
			vmax = np.nanmax(colors_valid)	

		# Create plot and add polygons to plot
		fig, ax = plt.subplots(figsize=(8, 8))
		norm = self.get_norm(data_yavg, norm_type, vmin=vmin, vmax=vmax)
		coll = PolyCollection(verts_valid, array=colors_valid, cmap=cmap,
			norm=norm, edgecolors=None)
		ax.add_collection(coll)

		# Plot LCFS
		if include_gfile:
			ax.contour(rgrid, zgrid, psirz, levels=[psi_lcfs], colors="k")
			ax.plot(rlim, zlim, color="k", lw=3)
			#ax.plot([2.1, 2.6], [Zdiv, Zdiv], linestyle="--", color="k")

		# Plot of the poloidal angle to help when referring where particles
		# are. Requires gfile to define where the magnetic axis is.
		if show_pol_ang and include_gfile:
			nang = 8
			line_len = 1.0
			for i in range(nang):
				angle = 2.0 * np.pi * i / nang  

				# End point of line
				R1 = R_axis + line_len * np.cos(angle)
				Z1 = Z_axis + line_len * np.sin(angle)

				ax.plot([R_axis, R1], [Z_axis, Z1], color="k", linestyle="--")

		# Add colorbar, show plot
		cbar = fig.colorbar(coll, ax=ax)
		cbar.set_label(data_name)
		ax.axis('equal')
		ax.autoscale_view()
		fig.show()
