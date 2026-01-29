import matplotlib.pyplot as plt
import netCDF4
import numpy as np
from matplotlib.colors import Normalize, LogNorm
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

		# Load data, correctly selecting the netCDF group it lives in
		# Geometry data are 3D vectors and don't have a time index
		if data_name in ["b_x", "b_y", "b_z", "jacobian", "gij_00", "gij_01",
			"gij_02", "gij_11", "gij_12", "gij_22"]:
			return self.nc["geometry"][data_name][:].data

		# Derived quantities
		elif data_name == "cyclotron_frequency":
			return self.calc_cyclo_freq()[frame]

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
		xlabel="x (m)", ylabel="y (m)", cbar_label=None, rsep=0.0,
		own_data=None):
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
		xlabel="x (m)", ylabel="z (m)", cbar_label=None, own_data=None):
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
		cbar_label=None, rsep=0.0, own_data=None):
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
		cbar_label=None, rsep=0.0, own_data=None):
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

	def calc_exb_drift(self):
		"""
		Calculate and return the X,Y,Z components of the ExB drift

		Outputs: Returns a tuple of the X,Y,Z ExB drift components
		"""
		
		# Pull out some arrays
		EX = self.nc["background"]["E_X"][:]
		EY = self.nc["background"]["E_Y"][:]
		EZ = self.nc["background"]["E_Z"][:]
		BX = self.nc["background"]["B_X"][:]
		BY = self.nc["background"]["B_Y"][:]
		BZ = self.nc["background"]["B_Z"][:]
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
	
	def calc_gradb_drift(self, charge):
		"""
		Calculate and return the X,Y,Z components of the grad-B drift

		Inputs
			- charge: Charge of the impurity

		Outputs: Returns a tuple of the X,Y,Z grad-B drift components
		"""
		
		# Pull out some arrays
		vX = self.nc["output"]["v_X"][:]
		vY = self.nc["output"]["v_Y"][:]
		vZ = self.nc["output"]["v_Z"][:]
		BX = self.nc["background"]["B_X"][:]
		BY = self.nc["background"]["B_Y"][:]
		BZ = self.nc["background"]["B_Z"][:]
		gradBX = self.nc["background"]["gradB_X"][:]
		gradBY = self.nc["background"]["gradB_Y"][:]
		gradBZ = self.nc["background"]["gradB_Z"][:]
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
	
	def calc_polarization_drift(self, charge):
		"""
		Calculate and return the X,Y,Z components of the polarization drift

		Inputs
			- charge: Charge of the impurity

		Outputs: Returns a tuple of the X,Y,Z polarization drift components
		"""

		# Pull out some arrays
		vX = self.nc["output"]["v_X"][:]
		vY = self.nc["output"]["v_Y"][:]
		vZ = self.nc["output"]["v_Z"][:]
		BX = self.nc["background"]["B_X"][:]
		BY = self.nc["background"]["B_Y"][:]
		BZ = self.nc["background"]["B_Z"][:]
		mz_kg = self.nc["input"]["imp_mass_amu"][:] * amu_to_kg
		t = self.nc["geometry"]["time"][:]

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
		return v_pol_X, v_pol_Y, v_pol_Z
	
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

