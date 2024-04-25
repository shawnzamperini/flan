# distutils: language=c++

# Copying the example found here:
# https://cython.readthedocs.io/en/stable/src/userguide/wrapping_CPlusPlus.html#wrapping-cplusplus
# Main particle following loop. Must be compiled with as cython code. 
# See setup.py file for details. 
import sys
import numpy as np
import constants
import time
import pickle
import random
from tqdm import tqdm
from scipy.interpolate import RegularGridInterpolator
import cython
from cython.parallel cimport prange

# Import Impurity C++ class.
from imp_cpp import PyImpurity

# OpenADAS interface (should convert to C++ eventually), but as of now
# it is normal python code that is compiled to C++. 
import openadas

from libcpp.vector cimport vector
from libcpp.string cimport string


# These seem to still be clocked as python objects when I use them, so
# I'll just define them at this scope instead.
cdef float ev = constants.ev

def follow_impurities(input_opts, gkyl_bkg, parallel=False):
    """
    This function contains the primary impurity ion following routine.
    """
    
    # Define variables to use as much C++ as possible here. Some more variables
    # are defined later since there is not strong convention to define them
    # all up front or as they're introduced. 
    cdef int i, j, imp_zstart_int, x_idx, y_idx, z_idx, fstart, fend, f
    cdef int loop_counts, iz_warn_count, rc_warn_count, coll_count
    cdef int var_count_qualifier, num_add_imps, coll_warn_count, charge
    cdef float x, y, z, local_ne, local_te, local_ti, local_ni
    cdef float local_ex, local_ey, local_ez, iz_coef, rc_coef, iz_prob
    cdef float rc_prob, ran, local_bz, print_percent, perc_done
    cdef float avg_time_followed, local_dtedx, local_dtedy, local_dtedz
    cdef float local_dtidx, local_dtidy, local_dtidz, local_viz
    cdef float local_stop_time, fstart_fl, fend_fl
    cdef bint coll_occured = 0
    
    # Extract input options for simulation as C++ types.
    cdef float imp_amu = input_opts["imp_mass"]
    cdef float imp_mass = input_opts["imp_mass"] * constants.amu
    cdef float imp_init_charge = input_opts["imp_init_charge"]
    cdef int num_imps = input_opts["num_imps"]
    cdef float imp_xmin = input_opts["imp_xmin"]
    cdef float imp_xmax = input_opts["imp_xmax"]
    cdef int gkyl_fstart = input_opts["gkyl_fstart"]
    cdef int gkyl_fend = input_opts["gkyl_fend"]
    cdef int imp_atom_num = input_opts["imp_atom_num"]
    cdef float imp_scaling_fact = input_opts["imp_scaling_fact"]
    cdef float imp_zstart_val = input_opts["imp_zstart_val"]
    cdef bint coll_forces = input_opts["coll_forces"]
    cdef bint collisions = input_opts["collisions"]
    cdef bint imp_xyz_tracks = input_opts["imp_xyz_tracks"]
    cdef bint var_reduction = input_opts["var_reduction"]
    cdef float var_reduction_freq = input_opts["var_reduction_freq"]
    cdef float var_reduction_min_weight = input_opts["var_reduction_min_weight"]
        
    # Check the both collisional forces and fully kinetic collisions are not
    # both on. In theory the latter captures the other. Eventually should
    # remove anything with coll_forces since it is a (probably bad) 
    # approximation.
    if coll_forces and collisions:
        print("Error! coll_forces and collisions cannot both be on. Choose one")
        sys.exit()

    # The gkeyll coordinates are the cell centers, but we want to know the
    # edges of each cells for our purposes.
    xwidth = (gkyl_bkg["x"][1] - gkyl_bkg["x"][0])
    x_bins = np.arange(gkyl_bkg["x"][0] - xwidth / 2, 
        gkyl_bkg["x"][-1] + xwidth / 2, xwidth)
    ywidth = (gkyl_bkg["y"][1] - gkyl_bkg["y"][0])
    y_bins = np.arange(gkyl_bkg["y"][0] - ywidth / 2, 
        gkyl_bkg["y"][-1] + ywidth / 2, ywidth)
    zwidth = (gkyl_bkg["z"][1] - gkyl_bkg["z"][0])
    z_bins = np.arange(gkyl_bkg["z"][0] - zwidth / 2, 
        gkyl_bkg["z"][-1] + zwidth / 2, zwidth)

    # Convert the imp_zstart option into an integer - just easier.
    imp_zstart_int = 0
    imp_zstart_opt = input_opts["imp_zstart_opt"]
    if imp_zstart_opt == "proportional_to_ni":
        imp_zstart_int = 1
        
        # Using the z distribution of the main ion density, create
        # CDFs that we will choose from later. This is tricky to
        # visualize in your head, but we are creating a normalized
        # PDF along the z direction for each x,y coordinate, and then
        # using cumsum to convert it to a CDF. Do it for each frame.
        zstart_cdf = np.zeros(gkyl_bkg["ni"].shape)
        for f in range(zstart_cdf.shape[0]):
            for i in range(zstart_cdf.shape[1]):
                for j in range(zstart_cdf.shape[2]):
                    zstart_cdf[f,i,j] = np.cumsum(gkyl_bkg["ni"][f,i,j] / 
                        gkyl_bkg["ni"][f,i,j].sum())
        
    elif imp_zstart_opt == "single_value":
        imp_zstart_int = 2
        
    else:
        print("Error: imp_zstart input not recognized: {}"
            .format(imp_zstart_opt))
        sys.exit()
    
    # Get the dimensions of the Gkeyll volume.
    gkyl_x = gkyl_bkg["x"] 
    gkyl_y = gkyl_bkg["y"]
    gkyl_z = gkyl_bkg["z"]
    cdef float gkyl_xmin = gkyl_bkg["x"].min()
    cdef float gkyl_xmax = gkyl_bkg["x"].max()
    cdef float gkyl_ymin = gkyl_bkg["y"].min()
    cdef float gkyl_ymax = gkyl_bkg["y"].max()
    cdef float gkyl_zmin = gkyl_bkg["z"].min()
    cdef float gkyl_zmax = gkyl_bkg["z"].max()
    
    # Assume constant time step between frames.
    cdef float dt = gkyl_bkg["time"][1] - gkyl_bkg["time"][0]
    print("Simulation timestep: {:.2e} s".format(dt))
    
    # Load the Gkeyll background into memory views, which have faster access. 
    cdef float[:,:,:,:] gkyl_ne = gkyl_bkg["ne"]
    cdef float[:,:,:,:] gkyl_te = gkyl_bkg["te"]
    cdef float[:,:,:,:] gkyl_ni = gkyl_bkg["ni"]
    cdef float[:,:,:,:] gkyl_ti = gkyl_bkg["ti"]
    cdef float[:,:,:,:] gkyl_elecx = gkyl_bkg["elecx"]
    cdef float[:,:,:,:] gkyl_elecy = gkyl_bkg["elecy"]
    cdef float[:,:,:,:] gkyl_elecz = gkyl_bkg["elecz"]
    cdef float[:,:,:,:] gkyl_viz = gkyl_bkg["viz"]
    cdef float[:,:,:] gkyl_b = gkyl_bkg["b"]
    
    # These need to be converted from eV/m to J/m. 
    cdef float[:,:,:,:] gkyl_dtedx = gkyl_bkg["dtedx"] * constants.ev
    cdef float[:,:,:,:] gkyl_dtedy = gkyl_bkg["dtedy"] * constants.ev
    cdef float[:,:,:,:] gkyl_dtedz = gkyl_bkg["dtedz"] * constants.ev
    cdef float[:,:,:,:] gkyl_dtidx = gkyl_bkg["dtidx"] * constants.ev
    cdef float[:,:,:,:] gkyl_dtidy = gkyl_bkg["dtidy"] * constants.ev
    cdef float[:,:,:,:] gkyl_dtidz = gkyl_bkg["dtidz"] * constants.ev
    
    # Can't cdef within a control statement like if, so define here and
    # potentially overwrite if the option is on.
    cdef float[:,:,:,:] stop_time = np.zeros(gkyl_bkg["ne"].shape, 
        dtype=np.float32)
    cdef alpha = 0.0
    cdef beta = 0.0
    
    if coll_forces:

        # Calculate the stopping time (Stangeby Eq. 6.35) ahead of time.
        # Calculate as python numpy array then assign to C array.
        ln_alpha = 15.0
        ion_amu = 2.04  # Hardcoding deuterium in.
        tmp_stop_time = 1.47e13 * imp_amu * gkyl_bkg["ti"] * \
            np.sqrt(gkyl_bkg["ti"] / ion_amu) / \
            ((1 + ion_amu / imp_amu) * gkyl_bkg["ni"] *  
            np.power(imp_atom_num, 2) * ln_alpha)
        stop_time = tmp_stop_time
        
        # Calculate the alpha and beta parameters that go in front of the 
        # electron and ion temperature gradient forces, respectively. 
        mu = imp_mass / (imp_mass + constants.mi)
        alpha = 0.71 * np.power(imp_atom_num, 2)
        beta = 3 * (mu + 5 * np.sqrt(2) * 
            np.square(imp_atom_num) * (1.1 * np.power(mu, 5/2) - 0.35 * 
            np.power(mu, 3/2)) - 1) / (2.6 - 2 * mu + 5.4 * np.square(mu))

        
    # I think it's probably more efficient to load all the random
    # numbers we'll need up front (where possible).
    cdef vector[float] xstart_rans
    cdef vector[float] ystart_rans
    cdef vector[float] zstart_rans_coarse
    cdef vector[float] zstart_rans_fine
    cdef vector[float] fstart_rans
    xstart_rans.reserve(num_imps)
    ystart_rans.reserve(num_imps)
    zstart_rans_coarse.reserve(num_imps)
    zstart_rans_fine.reserve(num_imps)
    fstart_rans.reserve(num_imps)
    xstart_rans = np.random.random(num_imps)
    ystart_rans = np.random.random(num_imps)
    zstart_rans_coarse = np.random.random(num_imps)
    zstart_rans_fine = np.random.random(num_imps)
    fstart_rans = np.random.random(num_imps)
    
    # Create OpenADAS object to load the rate coefficients. These are:
    #   oa: The OpenADAS class object. 
    #   rate_df_iz: The ionization rate coefficients
    #   rate_df_rc: The recombination rate coefficients 
    # The two rate_df are pandas DataFrames. We could use the oa class
    # to get the rate coefficients for a given ne, Te, charge, but I 
    # found it is way quicker to create RegularGridInterpolators for
    # each charge ahead of time that can be called for a given ne, Te.
    oa, rate_df_rc, rate_df_iz = load_adas(input_opts)
    
    def create_interps(df, atom_num):
        interps = []
        
        # Loop from 1 to atom_num+1 since the df is 1-indexed.
        for charge in range(1, atom_num+1):
            
            # Data is in log10, let's just make it normal units since
            # computers don't care. The units are converted to SI:
            #   Te: eV
            #   ne: cm-3 --> m-3
            #   rates: cm3/s --> m3/s
            tes = np.power(10, df.columns.values)
            nes = np.power(10, df.loc[charge].index.values) * 1e6 
            rates = np.power(10, df.loc[charge].values).T * 1e-6  

            # RegularGridInterpolater is recommended for our case.
            interp = RegularGridInterpolator((tes, nes), rates, 
                bounds_error=False, fill_value=None)
            interps.append(interp)
            
        return interps
    
    # Remember, these lists are 0-indexed, so for a specific charge number
    # we will want to index it with charge-1. The functions are called
    # with interps[charge-1](local_te, local_ne).
    interps_rc = create_interps(rate_df_rc, imp_atom_num)
    interps_iz = create_interps(rate_df_iz, imp_atom_num)
    
    if collisions:

        # Calculate the collision frequency arrays up front. We intentionally 
        # pass the non-cythonized arrays to make life easier in the function. 
        # Doesn't matter much since this is only called once and we use numpy.
        print("Calculating collision frequencies...")
        coll_freq_arr_py = calc_coll_freq_arr(gkyl_bkg["ne"], gkyl_bkg["ni"], 
            gkyl_bkg["te"], gkyl_bkg["ti"], input_opts["imp_mass"], imp_atom_num) 

    else:

        # Create array of zeros. Use the Gkeyll array shape and add a first 
        # dimension that is the size of the number of impurity charge states
        # since that's what shape array we get if collisions are on.
        coll_shape = list(gkyl_bkg["ne"].shape)
        coll_shape.insert(0, imp_atom_num)
        coll_freq_arr_py = np.zeros(coll_shape)

    # Assign nan's to zero. nan's occur I believe because Gkeyll arrays
    # have some extra zeros in them that results in nan's getting calculated.
    # There doesn't seem to be any effect on the simulation from them. 
    coll_freq_arr_py[np.isnan(coll_freq_arr_py)] = 0.0

    # Create memory view of it for faster access.
    cdef float[:,:,:,:,:] coll_freq_arr = np.array(coll_freq_arr_py, 
        dtype=np.float32)

    # The ending frame will always be the end of the simulation, this
    # is the first dimension in our arrays (0-indexed, not Gkeyll 
    # indexed).
    fend = gkyl_ne.shape[0]
    fend_fl = float(fend)
    
    # Arrays to keep track of the total particle weight passing through
    # each cell as well as the counts. This is the Monte Carlo way to
    # calculate density.
    imp_weight_arr = np.zeros(gkyl_bkg["ne"].shape)
    imp_count_arr = np.zeros(gkyl_bkg["ne"].shape)
    imp_vx_arr = np.zeros(gkyl_bkg["ne"].shape)
    imp_vy_arr = np.zeros(gkyl_bkg["ne"].shape)
    
    # Lists to record the impurity ion tracks. I am going to hardcode
    # a safegaurd in here for now, since I don't want to accidentally run
    # this and make a massive output file. I think you should be able to
    # accomplish any goals you need with 1,000 particles anyways.
    if imp_xyz_tracks:
        if num_imps > 1000:
            print("Warning! imp_xyz_tracks is set to True and over" +\
                " 1,000 particles were input. To avoid memory" +\
                " issues, num_imps is being reduced to 1,000. This" +\
                " is a hardcoded safeguard.")
            num_imps = 1000
        xyz_lists_x = [[] for i in range(num_imps)]
        xyz_lists_y = [[] for i in range(num_imps)]
        xyz_lists_z = [[] for i in range(num_imps)]

    # Determine after how often variance reduction should be recalculated in
    # terms of every X particles. This may seem like a verbose way to implement
    # this, but just modulo'ing the impurity number is not general enough. If 
    # you have 7 particles and var_reduction_freq = 0.1, then the statement:
    #   if i / num_imp == var_reduction_freq
    # will never execute. 
    cdef int var_red_part_freq = 0
    cdef int var_red_part_count = 0

    # Fringe and silly case, but nonetheless. 
    if num_imps * var_reduction_freq < 1:
        var_red_part_freq = 1

    # Round to nearest integer.
    else:
        var_red_part_freq = round(num_imps * var_reduction_freq)
   
    # Debugging dictionary containing some statistics.
    stat_dict = {"low_xbound_count":0, "high_xbound_count":0, 
        "low_zbound_count":0, "high_zbound_count":0, 
        "low_xbound_time":0.0, "coll_events":0, "total_imp_time":0.0}

    # Initialize to avoid compiler warning.
    local_dtedx = 0.0
    local_dtedy = 0.0
    local_dtedz = 0.0
    local_dtidx = 0.0
    local_dtidy = 0.0
    local_dtidz = 0.0

    # Initialize some loop variables.
    loop_counts = 0
    iz_warn_count = 0
    rc_warn_count = 0
    coll_warn_count = 0
    avg_time_followed = 0.0
    var_count_qualifier = 0
    num_add_imps = 0
    z = 0.0
    imp_foll_tstart = time.time()

    # Loop through tracking the full path of one impurity at a time.
    print("Beginning impurity following...")
    for i in tqdm(range(0, num_imps)):

        # Before starting an impurity, check if variance reduction needs to
        # be recalculated. We use a scheme based on splitting particles using
        # the collision probability in high variance (low count) regions. 
        if var_red_part_count == var_red_part_freq:

            # Use the median number of counts so far as the criteria for if a
            # particle is split (below median indicates the particle is in a 
            # low count region, which means it's important). 
            var_count_qualifier = np.median(
                imp_count_arr[np.nonzero(imp_count_arr)])
            print("var_count_qualifier: {}".format(var_count_qualifier))
            var_red_part_count = 0
        else:
            var_red_part_count += 1
 
        # Assume impurity can start at any frame (this means we are 
        # assuming a constant impurity source). This could lead to 
        # low statistics at the starting frames, should be studied. 
        fstart_fl = fend_fl * fstart_rans[i]
        fstart = np.floor(fstart_fl)

        # Determine random starting location. First is uniform between
        # xmin and xmax. 
        x = imp_xmin + (imp_xmax - imp_xmin) * xstart_rans[i]
        
        # We assume symmetry in the y direction, so uniformily
        # distributed between those as well.
        y = gkyl_ymin + (gkyl_ymax - gkyl_ymin) * ystart_rans[i]
        
        # Next determine starting z location. 
        #   1: This will be the closest value in the CDF, proportional
        #        to the ion density.
        #   2: Value specified. 
        if imp_zstart_int == 1:
            x_idx = np.argmin(np.abs(x - gkyl_x))
            y_idx = np.argmin(np.abs(y - gkyl_y))
            z_idx = np.argmin(np.abs(zstart_cdf[fstart, x_idx, y_idx] 
                - zstart_rans_coarse[i]))

            # To avoid picking at the same z location, chose uniformly
            # between each bounding z value.
            z = (gkyl_z[z_idx] - zwidth / 2) + zwidth * zstart_rans_fine[i]

        elif imp_zstart_int == 2:
            z = imp_zstart_val
                       
        # Now create our Impurity object. Starts with weight = 1.0.
        imp = PyImpurity(imp_atom_num, imp_mass, x, y, z, imp_init_charge, 
            fstart, 1.0)

        # Counter to track how many collisions the ion had.
        coll_count = 0

        # Put impurity into a list. At first this may seem confusing, why put 
        # it in a list by itself? But if particles are split then we would add 
        # that additional split particle to this list. The idea then is to
        # only move onto the next particle once this list is empty. You can
        # think of the number of particles as:
        # Particle #1
        # Particle #2 - Split
        #   Particle #2a - The split particle, added to list.
        # Particle #3 - Split
        #   Particle #3a - The split particle, added to list, which gets split 
        #                   again
        #       Particle #3b - The split particle, added to list
        # And so on...
        additional_imps = [imp]
        
        # Debug print statement to get a sense of where impurities are
        # being born.
        #if i > 0 and i % (num_imps / print_percent) == 0:
        #    print("  x = {:.4f}".format(x))
        #    print("  y = {:.4f}".format(y))
        #    print("  z = {:.4f}".format(z))
        #    print("  f = {:}".format(fstart))
        
        # Track transport of ions one frame at a time. Not all 
        # impurities will start at the first Gkeyll frame by design. The
        # first frame an impurity starts at is fstart, which when 
        # 0-indexed is fstart - gkyl_fstart, and the final frame it can
        # reach is gkyl_end-gkyl_fstart when 0-indexed.
        #for f in range(fstart, fend):
        while additional_imps:

            # Pop the ion from the list (this removes it from the list). 
            # Therefore the "while additional_imps:" loop will not execute
            # once this list is empty.
            imp = additional_imps.pop()
            for f in range(imp.fstart, fend):
                
                # Update nearest index values.
                x_idx = np.argmin(np.abs(imp.x - gkyl_x))
                y_idx = np.argmin(np.abs(imp.y - gkyl_y))
                z_idx = np.argmin(np.abs(imp.z - gkyl_z))

                # Add particle weight and count in our arrays. The particle
                # weight is multiplied by the timestep. You can think of this as
                # a particle in a cell is really only dt / (1 s) of a particle.
                #print(imp.weight)
                imp_weight_arr[f, x_idx, y_idx, z_idx] += imp.weight * dt
                imp_count_arr[f, x_idx, y_idx, z_idx] += 1
                
                # Likewise for particle velocity.
                imp_vx_arr[f, x_idx, y_idx, z_idx] += imp.vx
                imp_vy_arr[f, x_idx, y_idx, z_idx] += imp.vy
                
                # If desired, save the tracks of each impurity in (x,y,z)
                # and/or (vx,vy,vz) space.
                #if imp_xyz_tracks:
                #    xyz_lists_x[i].append(imp.x)
                #    xyz_lists_y[i].append(imp.y)
                #    xyz_lists_z[i].append(imp.z)
                
                # Extract the plasma parameters at the current location.
                local_ne = gkyl_ne[f, x_idx, y_idx, z_idx]
                local_te = gkyl_te[f, x_idx, y_idx, z_idx]
                local_ni = gkyl_ni[f, x_idx, y_idx, z_idx]
                local_ti = gkyl_ti[f, x_idx, y_idx, z_idx]
                local_ex = gkyl_elecx[f, x_idx, y_idx, z_idx]
                local_ey = gkyl_elecy[f, x_idx, y_idx, z_idx]
                local_ez = gkyl_elecz[f, x_idx, y_idx, z_idx]
                local_bz = gkyl_b[x_idx, y_idx, z_idx]
                local_viz = gkyl_viz[f, x_idx, y_idx, z_idx]
                local_stop_time = stop_time[f, x_idx, y_idx, z_idx]
                
                # Only need these if the collisional forces are on.
                if coll_forces:
                    local_dtedx = gkyl_dtedx[f, x_idx, y_idx, z_idx]
                    local_dtedy = gkyl_dtedy[f, x_idx, y_idx, z_idx]
                    local_dtedz = gkyl_dtedz[f, x_idx, y_idx, z_idx]
                    local_dtidx = gkyl_dtidx[f, x_idx, y_idx, z_idx]
                    local_dtidy = gkyl_dtidy[f, x_idx, y_idx, z_idx]
                    local_dtidz = gkyl_dtidz[f, x_idx, y_idx, z_idx]
                                    
                # Use ADAS to determine ionization/recombination 
                # probabilities, then pull a random number and see if we
                # are doing either of those. First load the rate 
                # coefficients by calling the function for this charge state
                # that we already loaded before the main loop.
                #print("te={:.2f}  ti={:.2f}  ne={:.2e}  charge={}".format(
                #    local_te, local_ti, local_ne, imp.charge))
                charge = imp.charge
                rc_coef = interps_rc[charge-1]((local_te, local_ne))
                iz_coef = interps_iz[charge-1]((local_te, local_ne))
                    
                # The probability of a ion ionizing/recombining during the 
                # timestep is then: 
                #   prob = coef [m3/s] * local_ne [m-3] * dt [s]
                iz_prob = iz_coef * local_ne * dt
                rc_prob = rc_coef * local_ne * dt
                
                # Pull random number, if it is less than the probability
                # then we have an event. Unfortunately we are pulling two
                # random numbers every loop here, so this could probably be
                # optimized a bit.
                ran = np.random.random()
                if ran < iz_prob and imp.charge < imp_atom_num:
                    imp.charge += 1
                #print("iz: ran={:.2e}  iz_prob={:.2e}".format(ran, iz_prob))
                ran = np.random.random()
                if ran < rc_prob and imp.charge > 0:
                    imp.charge -= 1
                #print("rc: ran={:.2e}  rc_prob={:.2e}".format(ran, rc_prob))
                
                # Check for wrong numbers?
                if imp.charge < 0:
                    print("Error: imp.charge < 0 ({})".format(imp.charge))
                if imp.charge > imp.imp_atom_num:
                    print("Error: imp.charge > {} ({})".format(imp.charge, 
                    imp.imp_atom_num))
                if iz_prob > 1:
                    #print("Error: iz_prob ({}) > 1".format(iz_prob))
                    iz_warn_count += 1
                if rc_prob > 1:
                    #print("Error: rc_prob ({}) > 1".format(rc_prob))
                    rc_warn_count += 1
                    
                # Perform step. imp object is updated within.
                #print("Before imp_step: {:.6f}  {:.6f}  {:.6f} ({:.3e} " \
                #   "{:.3e}  {:.3e})".format(imp.x, imp.y, imp.z, imp.vx, 
                #   imp.vy, imp.vz))
                imp_step(imp, local_ex, local_ey, local_ez, local_bz, dt, 
                    local_dtedx, local_dtedy, local_dtedz, local_dtidx, 
                    local_dtidy, local_dtidz, local_viz, local_stop_time, 
                    alpha, beta, coll_forces)
                #print("After imp_step: {:.6f}  {:.6f}  {:.6f} ({:.3e} " \
                #   "{:.3e}  {:.3e})".format(imp.x, imp.y, imp.z, imp.vx, 
                #   imp.vy, imp.vz))
                #print("{}-{}: x={:.4f}  y={:.4f}  z={:.4f}  ({}, {}, {})"
                #    .format(i, f, imp.x, imp.y, imp.z, x_idx, y_idx, z_idx))

                if collisions:
                    
                    # Test for a collision. If a collision occurs, the velocity
                    # vector will be shifted by 90 degrees in a random 
                    # direction.
                    ran = np.random.random()
                    charge = imp.charge
                    coll_freq = coll_freq_arr[charge-1, f, x_idx, y_idx, z_idx]
                    coll_prob = coll_freq * dt
                    #print("{:} {:} {:} {:} {:} Coll: {:.2e}  " \
                    #   "Prob = {:.3e}".format(imp.charge, f, x_idx, y_idx, 
                    #   z_idx, coll_freq, coll_prob))
                    if coll_prob > 1:
                        print("Error: coll_prob > 1 ({})".format(coll_prob))
                        coll_warn_count += 1

                    # Collision occured! Rotate velocity vector 90 degrees.
                    if ran < coll_prob:
                        rotate_vel_vector(imp)
                        coll_count += 1
                        stat_dict["coll_events"] += 1

                    # If a collision does not occur, we can still leverage the
                    # probability of a collision to split the particle in low
                    # count (high important) regions to improve statistics. 
                    # Sorry for the ugly if statements, I blame the deep
                    # indentation!
                    elif var_reduction and (imp.weight > 
                        var_reduction_min_weight):
                        if (imp_count_arr[f,x_idx,y_idx,z_idx] <=
                            var_count_qualifier):

                            # Split the particle. Each new particle weight is 
                            # proportional to the collision probability. One 
                            # will just be a continuation of the current 
                            # particle with reduced weight, while the other is 
                            # added to the list to be followed later. Its 
                            # weight is the original particle weight times the 
                            # collision probability. You can think of this as
                            # there is still only the original total weight to
                            # keep track of, now it's just split between two
                            # particles. The accounting all adds up to the same
                            # in the end.
                            imp.weight = imp.weight * (1 - coll_prob)
                            additional_imps.append(PyImpurity(imp_atom_num, 
                                imp_mass, imp.x, imp.y, imp.z, imp.charge, f, 
                                imp.weight * coll_prob))
                            num_add_imps += 1
              
                loop_counts += 1
                avg_time_followed += dt
                stat_dict["total_imp_time"] += dt
                
                # Bounds checking. Treat x and z bounds as absorbing, and
                # the y bound as periodic.
                if imp.x <= gkyl_xmin or imp.x >= gkyl_xmax:
                    if imp.x <= gkyl_xmin:
                        stat_dict["low_xbound_count"] += 1
                        stat_dict["low_xbound_time"] += imp.t
                    elif imp.x >= gkyl_xmax:
                        stat_dict["high_xbound_count"] += 1
                    break
                if imp.z <= gkyl_zmin or imp.z >= gkyl_zmax:
                    if imp.z <= gkyl_zmin:
                        stat_dict["low_zbound_count"] += 1
                    elif imp.z >= gkyl_zmax:
                        stat_dict["high_zbound_count"] += 1
                    break
                if imp.y <= gkyl_ymin:
                    imp.y = gkyl_ymax + (imp.y - gkyl_ymin)
                elif imp.y >= gkyl_ymax:
                    imp.y = gkyl_ymin + (imp.y - gkyl_ymax)
                
            # Free up memory (if necessary since it gets overwritten later).
            del imp

    # Calculate how long the main loop took.
    imp_foll_tend = time.time()
    cpu_time_used = imp_foll_tend - imp_foll_tstart
    print("Time spent following impurities: {:.1f} s".format(cpu_time_used))

    # Print how many additional impurities were printed.
    print("Additional split impurities followed: {}".format(num_add_imps))

    # Calculate the impurity density. This is only for 3D cartesian
    # grid. Normalize to the sum then divide by the volume of each cell
    # to get values in s/m3. Finally multiply by the scaling factor
    # (units of 1/s) to return the density in m-3.
    imp_dens_arr = imp_weight_arr / imp_weight_arr.sum()
    dx = gkyl_x[1] - gkyl_x[0]
    dy = gkyl_y[1] - gkyl_y[0]
    dz = gkyl_z[1] - gkyl_z[0]
    imp_dens_arr /= (dx * dy * dz)
    imp_dens_arr *= imp_scaling_fact
    
    # Divide by the counts in each cells for average velocity.
    with np.errstate(divide='ignore'):
        imp_vx_arr /= imp_count_arr
        imp_vy_arr /= imp_count_arr
    
    # Average time following an impurity ion.
    avg_time_followed /= (num_imps + num_add_imps)
    
    # Print how many time we had an ionization or recombination warning.
    if iz_warn_count > 0:
        print("Warning: The ionization probability (iz_prob) was " \
            ">1 {:} times, or {:.2e}% of the time. " \
            "Consider a smaller timestep.".format(iz_warn_count, 
            iz_warn_count/loop_counts*100))
    if rc_warn_count > 0:
        print("Warning: The recombination probability (rc_prob) was " \
            ">1 {:} times, or {:.2e}% of the time. " \
            "Consider a smaller timestep.".format(rc_warn_count, 
            rc_warn_count/loop_counts*100))
    if coll_warn_count > 0:
        print("Warning: The collision probability (coll_prob) was " \
            ">1 {:} times, or {:.2e}% of the time. " \
            "Consider a smaller timestep.".format(coll_warn_count, 
            coll_warn_count/loop_counts*100))

    # For now just print this, not gonna bother making it all pretty since
    # I'm abandoning this version of the code soon.
    print("stat_dict")
    print(stat_dict)

    # Bundle up the things we care about and return.
    return_dict = {
        "imp_dens_arr": imp_dens_arr, 
        "dt": dt, 
        "avg_time_followed": avg_time_followed, 
        "cpu_time_used": cpu_time_used,
        "imp_vx_arr": imp_vx_arr,
        "imp_vy_arr": imp_vy_arr,
        "stat_dict": stat_dict}
        
    return return_dict


# This decorator helps cut out the error checking done by python for 
# zero division. 
@cython.cdivision(True)
cdef imp_step(imp, float ex, float ey, float ez, float bz, float dt,
    float dtedx, float dtedy, float dtedz, float dtidx, float dtidy, 
    float dtidz, float viz, float stop_time, float alpha, float beta, 
    bint coll_forces):
    """
    Update impurity location. We are using the normal Lorentz force plus
    the friction, electron temperature gradient and ion temperature
    gradient forces in the z direction (these are the collisional
    forces). FF, FeG and FiG are defined in Stangeby, Eq. 6.21. 
    
    F = q(E + v x B) + Fc
    
    Fc = (see equation, too annoying to type out)
    """
    
    cdef float fx = 0.0
    cdef float fy = 0.0
    cdef float fz = 0.0
    cdef float dx = 0.0
    cdef float dy = 0.0
    cdef float dz = 0.0
    cdef float dvx = 0.0
    cdef float dvy = 0.0
    cdef float dvz = 0.0
    cdef float imp_charge = imp.charge
    cdef float imp_mass = imp.mass
    cdef float q = imp_charge * ev
    cdef float imp_vx = imp.vx
    cdef float imp_vy = imp.vy
    cdef float imp_vz = imp.vz
    
    # --------------------
    # Step in x direction.
    # --------------------
    # Electric field term
    fx += q * ex
    #print("fx_e {:.3e}".format(q * ex))
    
    # Magnetic field term
    fx += q * imp_vy * bz
    #print("fx_b {:.3e}".format(q * imp_vy * bz))
    
    # Collisional forces.
    if coll_forces:
        
        # Electron and ion temperature gradient forces.
        #print("x: feg = {:.3e}".format(alpha * dtedx))
        #print("x: fig = {:.3e}".format(beta * dtidx))
        #fx += alpha * dtedx
        #fx += beta * dtidx
        
        # Friction force.
        # Data not available from Gkeyll - need the cross field velocity 
        # moments.
    
        pass
    
    # --------------------
    # Step in y direction.
    # --------------------
    # Electric field term
    fy += q * ey
    #print("fy_e {:.3e}".format(q * ey))
    
    # Magnetic field term
    fy += -q * imp_vx * bz
    #print("fy_b {:.3e}".format(q * imp_vx * bz))
    
    if coll_forces:
        
        # Electron and ion temperature gradient terms.
        #print("y: feg = {:.3e}".format(alpha * dtedy))
        #print("y: fig = {:.3e}".format(beta * dtidy))
        #fy += alpha * dtedy
        #fy += beta * dtidy
        
        # Friction force.
        # Data not available from Gkeyll - need the cross field velocity 
        # moments.
        
        pass
    
    # --------------------
    # Step in z direction. 
    # --------------------
    # Electric field term
    fz += q * ez
    #print("fz_e {:.3e}".format(q * ez))
    
    if coll_forces:
        
        # Electron and ion temperature gradient terms.
        #print("z: feg = {:.3e}".format(alpha * dtedz))
        #print("z: fig = {:.3e}".format(beta * dtidz))
        fz += alpha * dtedz
        #print("fz_feg {:.3e}".format(alpha * dtedz))
        fz += beta * dtidz
        #print("fz_fig {:.3e}".format(beta * dtidz))
        
        # The friction force. 
        #print("z: ff = {:.3e}".format(imp_mass * (viz - imp_vz) / stop_time))
        #print("     viz = {:.3e}".format(viz))
        #print("     imp_vz = {:.3e}".format(imp_vz))
        #print("     stop_time = {:.3e}".format(stop_time))
        fz += imp_mass * (viz - imp_vz) / stop_time
        #print("fz_ff {:.3e}".format(imp_mass * (viz - imp_vz) / stop_time))
    
    # Now calculate the step sizes.
    dvx = fx * dt / imp_mass
    dvy = fy * dt / imp_mass
    dvz = fz * dt / imp_mass
    #print("  dvx = {:.2e}  dvy = {:.2e}  dvz = {:.2e}".format(dvx, 
    #    dvy, dvz))
    dx = (imp_vx + dvx) * dt
    dy = (imp_vy + dvy) * dt
    dz = (imp_vz + dvz) * dt
    #print("  dvx = {:.2e}  dvy = {:.2e}  dvz = {:.2e}".format(dvx, 
    #    dvy, dvz))
    
    
    # Update imp object with new values.
    imp.x += dx
    imp.y += dy
    imp.z += dz
    imp.vx += dvx
    imp.vy += dvy
    imp.vz += dvz
    imp.t += dt
    
    return None
    
    
def load_adas(input_opts):
    """
    Function to load the required OpenADAS files to keep track of 
    ionization/recombination on impurity ions. The files are:
    
    ACD: Effective recombination coefficients
    SCD: Effective ionization coefficients
    
    Inputs
    input_opts (dict): The input option dictionary.
    """
    
    # Select the correct OpenADAS file based on the impurity Z.  
    imp_atom_num = input_opts["imp_atom_num"]
    oa_dir = input_opts["openadas_dir"]
    if imp_atom_num == 74:
        acd_file = "acd50_w.dat"
        scd_file = "scd50_w.dat"
    else:
        print("Error: Files for imp_atom_num = {} not hardcoded " \
            "in yet!".format(imp_atom_num))
        sys.exit()
    acd_path = "{}/{}".format(oa_dir, acd_file)
    scd_path = "{}/{}".format(oa_dir, scd_file)
    
    
    # Load the OpenADAS data into pandas DataFrames.
    oa = openadas.OpenADAS()
    rate_df_rc = oa.read_rate_coef_unres(acd_path)
    rate_df_iz = oa.read_rate_coef_unres(scd_path)
    
    return oa, rate_df_rc, rate_df_iz

def calc_coll_freq_arr(ne, ni, te, ti, imp_mass, imp_atom_num):
    """
    Calculate and return an array of the impurity-ion collision frequency.
    Assumes impurity-deuterium collisions.
    
    Inputs
    ne : The electron density array from Gkeyll in m-3.
    ni : The main ion density array from Gkeyll in m-3.
    te : The electron temperature array from Gkeyll in eV.
    ti : The main ion temperature array from Gkeyll in eV.
    imp_mass : Impurity mass (amu).
    imp_atom_num : Impurity atomic number.

    Output
    A list where each element is a 4D array [f,x,y,z] of the collision
    frequency in order by charge state. So nu_imp_ion[0] is for the first
    impurity charge state, and so on.
    """
    
    # Load some constants up front. 
    amu = 1.660e-27  # amu to kg factor
    mi = 2.014  # amu
    mz = imp_mass # amu
    ev = 1.602e-19  # C
    eps0 = 8.85e-12
    Zi = 1  # Deuterium

    # Thermal velocities and reduced mass. v_fact is the factor out front such
    # that T is in eV and mass is in amu. 
    v_fact = np.sqrt(2 * ev / amu)
    vtz = v_fact * np.sqrt(ti / mz)
    vti = v_fact * np.sqrt(ti / mi)
    mr = mi * mz / (mi + mz)  # Still in amu.

    # This factor is all the constants out front of the lnalpha equation.
    lnalpha_fact = np.log(np.power(eps0, 3/2) * 4.0 * np.pi * amu / \
        np.power(ev, 5/2))

    # Likewise, this factor goes out front of the collision frequency equation.
    nu_corr_fact = np.power(ev, 4) * np.power(amu, 3/2) / (4 * np.pi * \
        np.square(eps0) * np.square(amu) * np.power(2 * ev, 3/2))

    # A list where each element is a 4D array [f,x,y,z] of the collision
    # frequency in order by charge state. So nu_imp_ion[0] is for the first
    # impurity charge state, and so on.
    nu_imp_ion = []

    # Loop across each atomic charge state (starting at 1). 
    for ZZ in tqdm(range(1, imp_atom_num+1)):
    
        # There is a small correction to lnalpha when considering impurity-ion
        # collisions. It's not major, e.g., for 10 eV and 1e19 m-3, it can be 
        # around 8. 
        #lnalpha_corr = np.log(np.sqrt(eps0 * te * ev / (ne * np.square(ev))) 
        #    * 4 * np.pi * eps0 * mr * np.square(vtz) / (ZZ * Zi 
        #    * np.square(ev)))
        lnalpha_corr = lnalpha_fact + np.log(np.sqrt(te / ne) * \
            np.square(vtz) * mr / (ZZ * Zi))

        # Collision frequency from Freidberg Eq. 9.48 with small correction
        # to account for charge of impurity ion and the aforementioned 
        # correction to lnalpha. 
        nu_corr = np.square(ZZ*ev) * np.square(Zi*ev) * ni * lnalpha_corr / \
            (4 * np.pi * np.square(eps0) * mz * mr) / (np.power(vtz, 3) + \
            1.3 * np.power(vti, 3))
        #num = np.square(ZZ*ev) * np.square(Zi*ev) * ni * lnalpha_corr
        #den1 = (4 * np.pi * np.square(eps0) * mz * mr)
        #den2 = (np.power(vtz, 3) + 1.3 * np.power(vti, 3))

        #print(nu_corr_fact)
        #print(ni[0,0,0,0])
        #print(lnalpha_corr[0,0,0,0])
        nu_corr = nu_corr_fact * np.square(ZZ) * np.square(Zi) * ni * \
            lnalpha_corr / (mz * mr) / (np.power(ti / mz, 3/2) + 1.3 * \
            np.power(ti / mi, 3/2))

        # Put into list.
        nu_imp_ion.append(nu_corr)

    # Need to return as numpy array for cython to easily convert to memoryview.
    return np.array(nu_imp_ion)

@cython.cdivision(True)
cdef rotate_vel_vector(imp):
    
    # Thankful for this stackoverflow post:
    # https://math.stackexchange.com/questions/1769501/how-do-i-rotate-a-vector-90-degrees-in-a-random-direction
    # Rotating by 90 degrees means the old and new vector will be perpendicular:
    #    a.b = x0x1 + y0y1 + z0z1 = 0
    # There is almost certainly a more efficient way to do this, but I'll
    # decribe a method this way.
    #  1. First randomly pick two of the components of the new vector (b).
    #  2. Assign those random values between 0-1. 
    #  3. The third component will then be, e.g.,:
    #       z1 = -(x0x1 + y0y1) / z0
    #  4. Normalize this new vector (b), and multiply by the length of the old
    #       vector. The result will be a vector that is perpindicular to our
    #       old vector. 

    # New vector to be filled in.
    cdef float[:] old_vect = np.array([imp.vx, imp.vy, imp.vz], 
        dtype=np.float32)
    cdef float[:] new_vect = np.array([0.0, 0.0, 0.0], dtype=np.float32)

    # Will need 4 random numbers. No dtype argument so making it double.
    cdef double [:] rans = np.random.random(2)

    # Pick two random direction and assign them random values. Do this by
    # shuffling 0, 1, 2 (corresponding to x, y, z), and then the first two get
    # random numbers and the last is the equation.
    cdef int[:] comps = np.array(random.sample([0, 1, 2], 3), dtype=np.intc)
    for i in range(2):
        new_vect[comps[i]] = rans[i]

    # Assign the third so it is perpendicular. Protect against divide by zero
    # errors by just assigning a small velocity.
    for i in range(len(old_vect)):
        if old_vect[i] == 0.0:
            old_vect[i] = 1e-30
    new_vect[comps[2]] = -(old_vect[comps[0]] * new_vect[comps[0]] + \
        old_vect[comps[1]] * new_vect[comps[1]]) / old_vect[comps[2]]

    # Normalize and multiply by magnitude of the original vector. We want to
    # chnage direction, not magnitude.
    cdef float old_vect_len = np.sqrt(old_vect[0]**2 + old_vect[1]**2 + \
        old_vect[2]**2)
    cdef float new_vect_len = np.sqrt(new_vect[0]**2 + new_vect[1]**2 + \
        new_vect[2]**2)

    # Assign the normalized components to the impurity's velocity vector.
    imp.vx = old_vect_len * new_vect[0] / new_vect_len
    imp.vy = old_vect_len * new_vect[1] / new_vect_len
    imp.vz = old_vect_len * new_vect[2] / new_vect_len

    #print("Old           New")
    #print("{:7.1f} {:7.1f}".format(old_vect[0], imp.vx))
    #print("{:7.1f} {:7.1f}".format(old_vect[1], imp.vy))
    #print("{:7.1f} {:7.1f}".format(old_vect[2], imp.vz))
    #dp = old_vect[0]*imp.vx + old_vect[1]*imp.vy + old_vect[2]*imp.vz
    #print("Dot product: {:.3e}".format(dp))




