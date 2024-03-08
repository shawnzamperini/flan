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
		
		# Load the plot coordinates, which are the same as Gkeyll's.
		x = self.nc["gkyl_bkg"]["gkyl_x"][:]
		y = self.nc["gkyl_bkg"]["gkyl_y"][:]
		z = self.nc["gkyl_bkg"]["gkyl_z"][:]
		
		# Find the closest index to the input z coordinate. 
		z_idx = np.argmin(np.abs(z - plot_z))
		
		# Load the impurity density at our z index.
		if dataname == "imp_dens":
			data = self.nc["imp_results"]["imp_dens"][:,:,:,z_idx]
		elif dataname == "ne":
			data = self.nc["gkyl_bkg"]["gkyl_ne"][:,:,:,z_idx]
		
		
		fig, ax = plt.subplots(figsize=(5, 4))
		
		# If a frame is specified, load that data, otherwise we'll be
		# making an .mp4 of all the plots.
		if f is not None:
			data_f = data[f].T

			# Create normalization based on input options.
			if normtype == "linear":
				norm = mpl.colors.Normalize(vmin=vmin, vmax=vmax)
			elif normtype == "log":
				norm = mpl.colors.LogNorm(vmin=vmin, vmax=vmax)

			# Make a 2D plot using pcolormesh. 
			cont = ax.pcolormesh(x, y, data_f, shading="nearest", 
				norm=norm, cmap=cmap)
			ax.set_xlabel("x (m)", fontsize=fontsize)
			ax.set_ylabel("y (m)", fontsize=fontsize)
			
			# Colorbar.
			div = make_axes_locatable(ax)
			cax = div.append_axes('right', '5%', '5%')
			cbar = fig.colorbar(cont, cax=cax)
			
			fig.tight_layout()
			fig.show()
		
			return data_f
		
		# Make an .mp4 of all the frames. Copying plot settings from 
		# above.
		else:
			data_f = data[0].T
	
			if normtype == "linear":
				norm = mpl.colors.Normalize(vmin=vmin, vmax=vmax)
			elif normtype == "log":
				norm = mpl.colors.LogNorm(vmin=vmin, vmax=vmax)
			cont = ax.pcolormesh(x, y, data_f, shading="nearest", 
				norm=norm, cmap=cmap)
			ax.set_xlabel("x (m)", fontsize=fontsize)
			ax.set_ylabel("y (m)", fontsize=fontsize)
			div = make_axes_locatable(ax)
			cax = div.append_axes('right', '5%', '5%')
			cbar = fig.colorbar(cont, cax=cax)
			cbar.set_label(r"$Impurity Density (\mathdefault{m^{-3}}$)", 
				fontsize=fontsize)
			fig.tight_layout()
			
			# Update function to pass to animation so it can create
			# our animation.
			def update(frame):
				
				# Clear axis before adding the new plot.
				ax.clear()
				
				# Chunk of code that is just a copy of above with the
				# comments removed. 
				data_f = data[frame].T
				if normtype == "linear":
					norm = mpl.colors.Normalize(vmin=vmin, vmax=vmax)
				elif normtype == "log":
					norm = mpl.colors.LogNorm(vmin=vmin, vmax=vmax)
				cont = ax.pcolormesh(x, y, data_f, shading="nearest", 
					norm=norm, cmap=cmap)
				ax.set_xlabel("x (m)", fontsize=fontsize)
				ax.set_ylabel("y (m)", fontsize=fontsize)
				# div = make_axes_locatable(ax)
				# cax = div.append_axes('right', '5%', '5%')
				
				# Clear colorbar axis.
				cax.cla()
				cbar = fig.colorbar(cont, cax=cax)
				cbar.set_label(r"Impurity Density ($\mathdefault{m^{-3}}$)", 
					fontsize=fontsize)
				fig.tight_layout()
	
				return cont
			
			# Create animation and save.
			if mp4name is None:
				mp4name = "{:}_imp_dens".format(self.nc["input_opts"]
				["case_name"][:][0])
			ani = animation.FuncAnimation(fig=fig, func=update, 
				frames=data.shape[0], interval=interval)
			writervideo = animation.FFMpegWriter(fps=15)
			ani.save('{:}.mp4'.format(mp4name), writer=writervideo)
			fig.show()
			
			
