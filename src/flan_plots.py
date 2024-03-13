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
		
	def plot_profiles(self, dataname, plot_z, plot_y=None, f=None, 
		fontsize=12, normtype="linear", vmin=None, vmax=None, 
		cmap="inferno", interval=80, fps=15, mp4name=None,
		linecolor="tab:red", plot_yscale="linear", plot_ymin=None,
		plot_ymax=None):
		"""
		Create 1D or 2D plot(s) of the impurity density at a given y 
		and/or z value.	If only plot_z is given, a 2D heatmap plot in 
		the xy plane is made. If plot_y is additionally given, then a 
		normal plot is made.
		
		Inputs
		dataname (str): 
		plot_z (float): z values at which the plot(s) are made. The 
			nearest z value will be chosen.
		plot_y (float): y values at which the plot(s) are made. The
			nearest y value will be chosen. If None (default) is given
			a 2D heatmap plot will be made. If a y value if given a
			standard plot will be made. A special keyword "average"
			will average along the y axis and also make a standard plot.
		f (int): Frame to make the plot at, 0-indexed. If None is 
			supplied, an .mp4 will be made instead of a static plot.
		fontsize (int): Fontsize for the plots.
		normtype (str): One of "linear" or "log". The scaling for the
			colormap.
		vmin (float): The minimum value on the colorbar.
		vmax (float): The maximum value on the colorbar.
		cmap (str): Colormap used in the plot. Uses matplotlib options.
		interval (int): Not sure what effect this has, may not matter 
			for mp4. 
		fps (int): The frames per second for the mp4. Smaller values 
			results in a slower playback speed. 
		mp4name (str): Name of the output mp4 video if you want it 
			different from the case name. Do not include .mp4 in the 
			name.
		linecolor (str): Matplotlib color for the line in the standard 
			plot when plot_y is specified.
		plot_yscale (str): Scale of the y-axis on the standard plot. 
			Valid options are "linear" or "log". If a list is specified,
			it must match the number of datanames where each entry 
			corresponds to a dataname entry.
		plot_ymin (float): The lower y-limit(s) of the standard 
			plot(s). Each entry corresponds to the dataname.
		plot_ymax (float): The upper y-limit(s) of the standard 
			plot(s). Each entry corresponds to the dataname.
		"""
		
		# Only acceptable string is average.
		if type(plot_y) is str and plot_y == "average":
			plot_y = -1
		elif type(plot_y) is str and plot_y != "average":
			print("Error: plot_y = {} not valid. ".format(plot_y) + \
			"Only valid string option is 'average'")
			return None
		
		# Conform to convention where this is a list, even if one entry.
		if type(dataname) is str:
			dataname = [dataname]
		if type(vmin) != list:
			vmin = [vmin]
		if type(vmax) != list:
			vmax = [vmax]
		if type(plot_yscale) != list:
			plot_yscale = [plot_yscale]
		if type(plot_ymin) != list:
			plot_ymin = [plot_ymin]
		if type(plot_ymax) != list:
			plot_ymax = [plot_ymax]	
		
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
		
		# Likewise for the y coordinate, if desired. Otherwise the data
		# will get averaged over the y coordinate. 
		if plot_y is None or str(plot_y) == "average":
			y_idx = np.arange(len(y), dtype=np.int)
		else:
			y_idx = np.argmin(np.abs(y - plot_y))
		
		# Load the data at our z index. Allow for multiple datanames at 
		# the same time for multiple plots. 
		data = []
		for dn in dataname:
			if dn == "nz":
				d = self.nc["imp_results"]["imp_dens"][:,:,y_idx,z_idx]
			elif dn == "ne":
				d = self.nc["gkyl_bkg"]["gkyl_ne"][:,:,y_idx,z_idx]
			elif dn == "te":
				d = self.nc["gkyl_bkg"]["gkyl_te"][:,:,y_idx,z_idx]
			elif dn == "ni":
				d = self.nc["gkyl_bkg"]["gkyl_ni"][:,:,y_idx,z_idx]
			elif dn == "ti":
				d = self.nc["gkyl_bkg"]["gkyl_ti"][:,:,y_idx,z_idx]
			elif dn == "ex":
				d = self.nc["gkyl_bkg"]["gkyl_elecx"][:,:,y_idx,z_idx]
			elif dn == "ey":
				d = self.nc["gkyl_bkg"]["gkyl_elecy"][:,:,y_idx,z_idx]
			elif dn == "ez":
				d = self.nc["gkyl_bkg"]["gkyl_elecz"][:,:,y_idx,z_idx]
		
			# If no y value input then average along the y axis 
			# (ignoring any zeros). 
			if str(plot_y) == "average":
				d[d == 0] = np.nan
				d = np.nanmean(d, axis=2)
				
			data.append(d)
		
		# Timestep used in the impuritty transport simulation.
		dt = float(self.nc["imp_results"]["imp_dt"][:])
		
		if len(data) == 1:
			fig, ax1 = plt.subplots(figsize=(5, 4))
			axs = [ax1]
		elif len(data) == 2:
			fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(9, 4))
			axs = [ax1, ax2]
		elif len(data) == 3:
			fig, (ax1, ax2, ax3) = plt.subplots(1, 2, figsize=(11, 4))
			axs = [ax1, ax2, ax3]
		elif len(data) == 4:
			fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(8, 8))
			axs = [ax1, ax2, ax3, ax4]
		
		if plot_y is None:
			caxs = []
			for ax in axs:
				
				# A grey background looks better generally.
				ax.set_facecolor("grey")
		
				# Colorbars.
				div = make_axes_locatable(ax)
				caxs.append(div.append_axes('right', '5%', '5%'))
		
		# If a frame is specified, load that data, otherwise we'll be
		# making an .mp4 of all the plots.
		if f is not None:
			
			# Assemble one plot as a time.
			for i in range(0, len(data)):
				
				# 2D heatmap.
				if plot_y is None:
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
				
				# Standard plot.
				else:
					data_f = data[i][f]
					axs[i].plot(x, data_f, color=linecolor, lw=2)
					axs[i].set_xlabel("x (m)", fontsize=fontsize)
					axs[i].set_title("t = {:.2f} us".format(dt * f * 1e6))
					axs[i].set_yscale(plot_yscale[i])
					axs[i].set_ylim(plot_ymin[i], plot_ymax[i])
			
			fig.tight_layout()
			fig.show()
		
			return data_f
		
		# Make an .mp4 of all the frames. Copying plot settings from 
		# above.
		else:
			
			# Assemble one plot at a time.
			for i in range(0, len(data)):
				
				# 2D heatmap.
				if plot_y is None:
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
				
				# Standard plot.
				else:
					data_f = data[i][0]
					axs[i].plot(x, data_f, color=linecolor, lw=2)
					axs[i].set_xlabel("x (m)", fontsize=fontsize)
					axs[i].set_title("t = {:.2f} us".format(0))
					axs[i].set_yscale(plot_yscale[i])
					axs[i].set_ylim(plot_ymin[i], plot_ymax[i])
				
				fig.tight_layout()
			
			# Update function to pass to animation so it can create
			# our animation.
			if plot_y is None:
				
				# The heatmap update.
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
			
			else:
				
				# The standard plot update.
				def update(frame):
					data_f = data[i][frame]
					axs[i].clear()
					axs[i].plot(x, data_f, color=linecolor, lw=2)
					axs[i].set_xlabel("x (m)", fontsize=fontsize)
					axs[i].set_title("t = {:.2f} us".format(dt * frame * 1e6))
					axs[i].set_yscale(plot_yscale[i])
					axs[i].set_ylim(plot_ymin[i], plot_ymax[i])
			
			# Create animation and save.
			if mp4name is None:
				mp4name = "{:}_imp_dens".format(self.nc["input_opts"]
				["case_name"][:][0])
			ani = animation.FuncAnimation(fig=fig, func=update, 
				frames=data[0].shape[0], interval=interval)
			writervideo = animation.FFMpegWriter(fps=fps)
			ani.save('{:}.mp4'.format(mp4name), writer=writervideo)
			fig.show()
			
			
