# To be included in flan_plots.py potentially? Will need to access the data from flan_plots



import numpy as np
import matplotlib as plt









# Takes a frame and radial index. Calculates the average velocity from non-zero velocity particles in a given frame at that radial index. Note that there are 96 radial indices. 

def frame_avg_velocity(radialIndex, frame_index):

    # Gets the impurity velocity (imp_vX) data from flan_plots.py
    a, b, c = fp.plot_frame_xy("imp_vX", frame=frame_index, z0=0.01, showplot=False, norm_type="log")

    # Initialize the variables to be used. Total is just the total counter for the velocity, and velocity_average is the velocity average that will be calculated after going through the length of the data array
    
    total = 0.0
    count = 0.0
    velocity_average = 0.0

    # This loop will go through the length of the velocity array, and add each non-zero velocity to total. The counter will also increase by one for each iterated non-zero velocity, 
    # thus keeping track of how many velocities we are averaging over. 

    for i in range(len(c)):  # i from 0 to 95
        if(c[i][radialIndex] != 0) 
            total += c[i][radialIndex]
            count +=1

    # Now, calculate the average velocity and return it! Returns zero in the case that the radial span has no particles - thus giving a count=0
    if count > 0:
        velocity_average = total / count
    else:
        velocity_average = 0.0

    return velocity_average


# This function takes an array of 2D coordinates and 1D of data, then gives the poloidally averaged value of the data at a single radial index. The radial index must likewise be specified.

def radial_index_avg_density(radialIndex, arr):
    total_radial_density=0.0
    average_radial_density=0
    count=0

    # Wanted to define a function to find the average density at a given radial coordinate, since we need to use it for time step i and i+1 to calculate a derivative. This loops through the poloidal coordinate, 
    # summing to get a total density for a given radial coordinate. Then, it returns the average by dividing by the total length of the array. 
        
    for i in range(arr.shape[0]):
        
        # Checking to make sure the entry is not zero. We only calculate averages for non-zero densities, thus introducing the non-zero density counter. 
                
        if arr[radialIndex][i][8] != 0: 
            count+=1
            total_radial_density+=arr[radialIndex][i][8]
                
    # Returns the average, if the count is greater than one. It should be... but it doesn't hurt to add the failsafe. 

    if count > 0:
        average_radial_density = total_radial_density / count
    else:
        average_radial_density = 0.0

    
    return average_radial_density


# Takes the change in density from the current radial index to the next one, and divides it by the spatial difference between radial indices. 

def frame_density_derivative(frame_index, radial_index):

    # Inputs: Takes in the index of the frame you'd like to look at, and takes the index of the radial coordinate. 
    # Outputs: Gives you the derivative of the density between your radial index and the next radial index 


    # Begin by getting the 3D array of 2D spatial coordinates and density as the 3rd. 
    density_array = fp.load_data_frame("imp_density",frame=frame_index)

    
    # Finds the poloidal averages of density at the given radial index and that at the radial index +1. Then we find the derivative. When this is called and you're looking at the last frame in the simulation,
    # instead of throwing an index error, you'll just get zero. 
    
    avg_now = radial_index_avg_density(radial_index, density_array) 

    if radial_index != density_array.shape[0]:
        avg_next = radial_index_avg_density(radial_index + 1, density_array) 
    else: 
        avg_next=avg_now

    # The difference between the radial coordinates is held around 0.0005746 within some small machine precision amount. 
    
    derivative = (avg_next-avg_now)/0.0005746

    return derivative


# We now have everything necessary to compute the diffusion coefficient - so I'll go ahead and do that. 
def frame_radial_diffusion(radial_index, frame_index): 
    
    # Want to get the density array once more, since calculation of the diffusion coefficient involves a local density. 
    
    density_array = fp.load_data_frame("imp_density",frame=frame_index)
    
    #define diffusion coefficient. 

    diffusion_coefficient = 0.0

    # Now, we can calculate each component in terms of other functions used. 

    numerator_diffusion_coefficient = -frame_avg_velocity(radialIndex=radial_index, frame_index) * radial_index_avg_density(radial_index, density_array)

    denominator_diffusion_coefficient = frame_density_derivative(frame_index, radial_index)

    # Last, we calculate the coefficient and return it. The density derivative is zero at the boundary point because of how it was discretized, so 
    # the if-statement avoids the divide-by-zero error by just setting the diffusivity equal to 0. I'm very open to other solutions to mitigate this. 

    if denominator_diffusion_coefficient != 0:
        diffusion_coefficient = numerator_diffusion_coefficient/denominator_diffusion_coefficient
    else: 
        diffusion_coefficient = 0

    
    return diffusion_coefficient 




# Now we want to create an array of diffusion coefficients. There will be one diffusion coefficient for each radial index. This array will be used in the plots. 

def frame_diffusion_array(v_frame_index)

    # This will be the array that's returned. It's as long as the radial component array. 

    diffusion_coef_array=np.zeros((fp.load_data_frame("imp_density",frame=v_frame_index)).shape[0])

    # Now we iterate over all radial lengths, calculate the diffusion coefficient, and put them in the array
    
    for i in range((fp.load_data_frame("imp_density",frame=v_frame_index)).shape[0])

        # We want to fill an array with the values of the diffusion coefficient based on the radial length, so we let the iterated variable be the radial length
        
        diffusion_coef_array[i]=frame_radial_diffusion(radial_index=i, frame_index = v_frame_index)

    return diffusion_coef_array





# Finally, we want to plot the  diffusion coefficient vs. the radial coordinates. Questionable if this function should be in this file or a parent file, since we'll 
# eventually want to plot multiple species against each other and this framework won't allow for that. 

def radial_diffusion_plot(i_frame_index, species) 

    # We need to get the array of x-values now, and then put them in the plot. 

    x,y,data=fp.plot_frame_xy("imp_density",frame=199,z0=0.01,norm_type="log",vmin=1e-5,vmax=1e-1)
    plt.plot(x, frame_diffusion_array(i_frame_index))

    # Now we'll just label the plot's axes and specify the species. 
    
    plt.xlabel("Radial length (x)")
    plt.ylabel("Diffusion Coefficient (y)")
    plt.title("Radial Diffusion for " + species)

    # Show the plot!
    
    plt.grid(True)
    plt.show()





























