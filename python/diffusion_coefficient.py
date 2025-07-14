import numpy as np
import flan_plots
fp=flan_plots.FlanPlots(self.nc)

delta_x = 0.0005746

# Takes a frame and radial index. Calculates the average velocity from non-zero entries in a given frame at that radial index. Note that there are 96 radial indices. 
def frame_avg_velocity(radialIndex, frame_index):
    a, b, c = fp.plot_frame_xy("imp_vX", frame=frame_index, z0=0.01, showplot=False, norm_type="log")
    
    total = 0.0
    for i in range(96):  # i from 0 to 95
        if(c[i][radialIndex != 0) 
            total += c[i][radialIndex]
    
    velocity_average = total / 96.0
    return velocity_average






# Takes the change in density from the current radial index to the next one, and divides it by the spatial difference between radial indices. 
def frame_density_derivative(frame_index, radial_index):

    # Inputs: Takes in the index of the frame you'd like to look at, and takes the index of the radial coordinate. 
    # Outputs: Gives you the derivative of the density between your radial index and the next radial index 


    # Begin by getting the 3D array of 2D spatial coordinates and density as the 3rd. 
    density_array = fp.load_data_frame("imp_density",frame=frame_index)

    def radial_index_avg_velocity(radialIndex):
        total_radial_density=0.0
        count=0

        # Wanted to define a function to find the average density at a given radial coordinate, since we need to use it for time step i and i+1 to calculate a derivative. This loops through the poloidal coordinate, 
        # summing to get a total density for a given radial coordinate. Then, it returns the average by dividing by the total length of the array. 
        
        for i in range(density_array.shape[1]):
        
            # Checking to make sure the entry is not zero. We only calculate averages for non-zero densities, thus introducing the non-zero density counter. 
                
            if density_array[radialIndex][i][8] != 0: 

                
                count++
                total_radial_density+=density_array[radialIndex][i][8]

        
        return total_radial_density/count
                
        
    # Finds the poloidal averages at the given radial index and that at the radial index +1. Then we find the derivative. 
    avg_now = radial_index_avg_velocity(radial_index) 
    avg_next = radial_index_avg_velocity(radial_index + 1) 
        
    derivative = (avg_next-avg_now)/0.0005746

    return derivative


# We now have everything necessary to compute the diffusion coefficient - so I'll go ahead and do that. 
def frame_diffusion(radial_index, frame): 




