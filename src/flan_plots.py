# Class object to interface with Flan output and help make plots.
from netCDF4 import Dataset
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from mpl_toolkits.axes_grid1 import make_axes_locatable
import numpy as np

class FlanPlots:
	
	def __init__(self, path):
		"""
		Inputs
		path (str): Path to .nc file from a Flan run.
		"""
		
		# Load the netCDF file and some of the simulation details.
		self.nc = Dataset(path)
		
	def __repr__(self):
		
		repr_str = "{}\n{}".format(self.nc.description, 
			self.nc.history)
		return repr_str
		
	def plot_xy(self, dataname, plot_z, f=None, fontsize=12, 
		normtype="linear", vmin=None, vmax=None, cmap="inferno",
		interval=80, mp4name=None):
		"""
		Create 2D plot(s) of the impurity density at a given z value.
		The plots will be in the xy plane.
		
		Inputs
		z (float): z values at which the plot(s) are made.
		f (int): Frame to make the plot at, 0-indexed. If None is 
			supplied, an .mp4 will be made instead of a static plot.
		"""
		
		# Conform to convention where this is a list, even if one entry.
		if type(dataname) is str:
			dataname = [dataname]
		if type(vmin) != list:
			vmin = [vmin]
		if type(vmax) != list:
			vmax = [vmax]
		
		if len(vmin) != len(dataname):
			print("Error! Please enter as many vmins as there are " \
				"datanames, as a list.")
			return None
		if len(vmax) != len(dataname):
			print("Error! Please enter as many vmaxs as there are " \
				"datanames, as a list.")
			return None
		
		# Only allow up to 4 plots.
		if len(dataname) > 4:
			print("Slow down there cowboy! Only 4 plots allowed " \
				"at a time.")
			return None
		
		# Load the plot coordinates, which are the same as Gkeyll's.
		x = self.nc["gkyl_bkg"]["gkyl_x"][:]
		y = self.nc["gkyl_bkg"]["gkyl_y"][:]
		z = self.nc["gkyl_bkg"]["gkyl_z"][:]
		
		# Find the closest index to the input z coordinate. 
		z_idx = np.argmin(np.abs(z - plot_z))
		
		# Load the impurity density at our z index. Allow for multiple
		# datanames at the same time for multiple plots. 
		data = []
		for dn in dataname:
			if dn == "imp_dens":
				data.append(self.nc["imp_results"]["imp_dens"]
					[:,:,:,z_idx])
			elif dn == "ne":
				data.append(self.nc["gkyl_bkg"]["gkyl_ne"][:,:,:,z_idx])
			elif dn == "te":
				data.append(self.nc["gkyl_bkg"]["gkyl_te"][:,:,:,z_idx])
			elif dn == "ni":
				data.append(self.nc["gkyl_bkg"]["gkyl_ni"][:,:,:,z_idx])
			elif dn == "ti":
				data.append(self.nc["gkyl_bkg"]["gkyl_ti"][:,:,:,z_idx])
			elif dn == "ex":
				data.append(self.nc["gkyl_bkg"]["gkyl_elecx"]
					[:,:,:,z_idx])
			elif dn == "ey":
				data.append(self.nc["gkyl_bkg"]["gkyl_elecy"]
					[:,:,:,z_idx])
			elif dn == "ez":
				data.append(self.nc["gkyl_bkg"]["gkyl_elecz"]
					[:,:,:,z_idx])
		
		# Timestep used in the impuritty transport simulation.
		dt = float(self.nc["imp_results"]["imp_dt"][:])
		
		if len(data) == 1:
			fig, ax1 = plt.subplots(figsize=(5, 4))
			axs = [ax1]
		elif len(data) == 2:
			fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(9, 4))
			ax2.set_facecolor("grey")
			axs = [ax1, ax2]
		elif len(data) == 3:
			fig, (ax1, ax2, ax3) = plt.subplots(1, 2, figsize=(11, 4))
			ax2.set_facecolor("grey")
			ax3.set_facecolor("grey")
			axs = [ax1, ax2, ax3]
		elif len(data) == 4:
			fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(8, 8))
			ax2.set_facecolor("grey")
			ax3.set_facecolor("grey")
			ax4.set_facecolor("grey")
			axs = [ax1, ax2, ax3, ax4]
		ax1.set_facecolor("grey")
		
		# Colorbars.
		caxs = []
		for i in range(0, len(data)):
			div = make_axes_locatable(axs[i])
			caxs.append(div.append_axes('right', '5%', '5%'))
		
		# If a frame is specified, load that data, otherwise we'll be
		# making an .mp4 of all the plots.
		if f is not None:
			
			# Assemble one plot as a time.
			for i in range(0, len(data)):
				data_f = data[i][f].T

				# Create normalization based on input options.
				if normtype == "linear":
					norm = mpl.colors.Normalize(vmin=vmin[i], 
						vmax=vmax[i])
				elif normtype == "log":
					norm = mpl.colors.LogNorm(vmin=vmin[i], 
						vmax=vmax[i])

				# Make a 2D plot using pcolormesh. 
				cont = axs[i].pcolormesh(x, y, data_f, shading="nearest", 
					norm=norm, cmap=cmap)
				axs[i].set_xlabel("x (m)", fontsize=fontsize)
				axs[i].set_ylabel("y (m)", fontsize=fontsize)
				axs[i].set_title("t = {:.2f} us".format(dt * f * 1e6))
				
				# Colorbar.
				#div = make_axes_locatable(axs[i])
				#cax = div.append_axes('right', '5%', '5%')
				cbar = fig.colorbar(cont, cax=caxs[i])
			
			fig.tight_layout()
			fig.show()
		
			return data_f
		
		# Make an .mp4 of all the frames. Copying plot settings from 
		# above.
		else:
			# Assemble one plot as a time.
			for i in range(0, len(data)):
				data_f = data[i][0].T

				# Create normalization based on input options.
				if normtype == "linear":
					norm = mpl.colors.Normalize(vmin=vmin[i], 
						vmax=vmax[i])
				elif normtype == "log":
					norm = mpl.colors.LogNorm(vmin=vmin[i], 
						vmax=vmax[i])

				# Make a 2D plot using pcolormesh. 
				cont = axs[i].pcolormesh(x, y, data_f, shading="nearest", 
					norm=norm, cmap=cmap)
				axs[i].set_xlabel("x (m)", fontsize=fontsize)
				axs[i].set_ylabel("y (m)", fontsize=fontsize)
				axs[i].set_title("t = {:.2f} us".format(0))
				
				# Colorbar.
				#div = make_axes_locatable(axs[i])
				#cax = div.append_axes('right', '5%', '5%')
				cbar = fig.colorbar(cont, cax=caxs[i])
			
				fig.tight_layout()
			
			# Update function to pass to animation so it can create
			# our animation.
			def update(frame):
				
				for i in range(0, len(data)):
					data_f = data[i][frame].T
					axs[i].clear()
					caxs[i].clear()

					# Create normalization based on input options.
					if normtype == "linear":
						norm = mpl.colors.Normalize(vmin=vmin[i], 
							vmax=vmax[i])
					elif normtype == "log":
						norm = mpl.colors.LogNorm(vmin=vmin[i], 
							vmax=vmax[i])

					# Make a 2D plot using pcolormesh. 
					cont = axs[i].pcolormesh(x, y, data_f, shading="nearest", 
						norm=norm, cmap=cmap)
					axs[i].set_xlabel("x (m)", fontsize=fontsize)
					axs[i].set_ylabel("y (m)", fontsize=fontsize)
					axs[i].set_title("t = {:.2f} us".format(dt * frame * 1e6))
					
					# Colorbar.
					#div = make_axes_locatable(axs[i])
					#cax = div.append_axes('right', '5%', '5%')
					cbar = fig.colorbar(cont, cax=caxs[i])
				
					fig.tight_layout()
				return cont
			
			# Create animation and save.
			if mp4name is None:
				mp4name = "{:}_imp_dens".format(self.nc["input_opts"]
				["case_name"][:][0])
			ani = animation.FuncAnimation(fig=fig, func=update, 
				frames=data[0].shape[0], interval=interval)
			writervideo = animation.FFMpegWriter(fps=15)
			ani.save('{:}.mp4'.format(mp4name), writer=writervideo)
			fig.show()
			
			
	def plot_x(self):
		"""
		Plot a quantity along the x coordinate. 
		"""
		pass
