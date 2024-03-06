# distutils: language=c++

# Copying the example found here:
# https://cython.readthedocs.io/en/stable/src/userguide/wrapping_CPlusPlus.html#wrapping-cplusplus
# Main particle following loop. Must be compiled with as cython code. 
# See setup.py file for details. 
import sys
import numpy as np
import constants

# Import Impurity C++ class.
from imp_cpp import PyImpurity

# OpenADAS interface (should convert to C++ eventually), but as of now
# it is normal python code that is compiled to C++. 
import openadas

from libcpp.vector cimport vector


def follow_impurities(input_opts, gkyl_bkg):
    """
    This function contains the primary impurity ion following routine.
    """
    
    # Define variables to use as much C++ as possible here.
    cdef int i, j, imp_zstart_int, x_idx, y_idx, z_idx, fstart
    cdef float x, y, z, local_ne, local_te, local_ti, local_ni
    cdef float local_ex, local_ey, local_ez, iz_coef, rc_coef, iz_prob
    cdef float rc_prob, ran, local_bz
    
    # Extract input options for simulation as C++ types.
    cdef float imp_mass = input_opts["imp_mass"] * constants.amu
    cdef float imp_init_charge = input_opts["imp_init_charge"]
    cdef int num_imps = input_opts["num_imps"]
    cdef float imp_xmin = input_opts["imp_xmin"]
    cdef float imp_xmax = input_opts["imp_xmax"]
    cdef int gkyl_fstart = input_opts["gkyl_fstart"]
    cdef int imp_atom_num = input_opts["imp_atom_num"]
    
    # If we interpolated frames, gkyl_fend will be that many additional
    # frames longer! 
    # To-do...
    cdef int gkyl_fend = input_opts["gkyl_fend"]
    
    # Other options where it isn't a big deal if they aren't C types.
    imp_zstart_opt = input_opts["imp_zstart_opt"]
    
    # Convert the imp_zstart option into an integer - just easier.
    if imp_zstart_opt == "proportional_to_ni":
        imp_zstart_int = 1
        
        # Using the z distribution of the main ion density, create
        # CDFs that we will choose from later. This is tricky to
        # visualize in your head, but we are creating a normalized
        # PDF along the z direction for each x,y coordinate, and then
        # using cumsum to convert it to a CDF. 
        zstart_cdf = np.zeros(gkyl_bkg["ni_avg"].shape)
        for i in range(zstart_cdf.shape[0]):
            for j in range(zstart_cdf.shape[1]):
                zstart_cdf[i,j] = np.cumsum(gkyl_bkg["ni_avg"][i,j] / 
                    gkyl_bkg["ni_avg"][i,j].sum())
        
    else:
        print("Error: imp_zstart input not recognized: {}"
            .format(imp_zstart_opt))
        sys.exit()
    
    # Get the dimensions of the Gkeyll volume.
    #nx = len(gkyl_bkg["x"])
    #ny = len(gkyl_bkg["y"])
    #nz = len(gkyl_bkg["z"])
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
    
    # Load the Gkeyll background into vectors. 
    # I'm not sure I fully grasp this... see here:
    # https://cython.readthedocs.io/en/stable/src/userguide/memoryviews.html
    # For now stick with pyhton numpy arrays but this is a potential
    # area for speed up.
    
    # I think it's probably more efficient to load all the random
    # numbers we'll need up front (where possible).
    cdef vector[float] xstart_rans
    cdef vector[float] ystart_rans
    cdef vector[float] zstart_rans
    cdef vector[float] fstart_rans
    xstart_rans.reserve(num_imps)
    ystart_rans.reserve(num_imps)
    zstart_rans.reserve(num_imps)
    fstart_rans.reserve(num_imps)
    xstart_rans = np.random.random(num_imps)
    ystart_rans = np.random.random(num_imps)
    zstart_rans = np.random.random(num_imps)
    fstart_rans = np.random.random(num_imps)
    
    # Create OpenADAS object to load the rate coefficients. These are:
    #   oa: The OpenADAS class object. 
    #   rate_df_iz: The ionization rate coefficients
    #   rate_df_rc: The recombination rate coefficients 
    # The two rate_df are pandas DataFrames. The oa class facilitates
    # extracting rate coefficients from the DataFrames via, e.g.:
    #   oa.get_rate_coef(rate_df_iz, te=15, ne=5e18, charge=5)
    # It will interpolate for values of te/ne that are not in the
    # DataFrame. The returned units are m3/s.
    oa, rate_df_rc, rate_df_iz = load_adas(input_opts)
    
    # Will need these available as floats.
    cdef float gkyl_fstart_fl = float(input_opts["gkyl_fstart"])
    cdef float gkyl_fend_fl = float(input_opts["gkyl_fend"])
    cdef float fstart_fl
    
    # Arrays to keep track of the total particle weight passing through
    # each cell as well as the counts. This is the Monte Carlo way to
    # calculate density.
    imp_weight_arr = np.zeros(gkyl_bkg["b"].shape)
    imp_count_arr = np.zeros(gkyl_bkg["b"].shape)
    
    # Printout every X percent.
    cdef float print_percent = 10

    # Loop through tracking the full path of one impurity at a time.
    for i in range(0, num_imps):
        
        # Print out counter every 10%. 
        if i > 0 and i % (num_imps / print_percent) == 0:
            print("{}/{}".format(i, num_imps))
        
        # Determine random starting location. First is uniform between
        # xmin and xmax. 
        x = imp_xmin + (imp_xmax - imp_xmin) * xstart_rans[i]
        x_idx = np.argmin(np.abs(x - gkyl_x))
        
        # We assume symmetry in the y direction, so uniformily
        # distributed between those as well.
        y = gkyl_ymin + (gkyl_ymax - gkyl_ymin) * ystart_rans[i]
        y_idx = np.argmin(np.abs(y - gkyl_y))
        
        # Next determine starting z location. This will be the closest
        # value in the CDF.
        if imp_zstart_int == 1:
            z_idx = np.argmin(np.abs(zstart_cdf[x_idx, y_idx] 
            - zstart_rans[i]))
        z = gkyl_z[z_idx]
        
        # Assume impurity can start at any frame (this means we are 
        # assuming a constant impurity source). This could lead to 
        # low statistics at the starting frames, should be studied. 
        fstart_fl = gkyl_fstart_fl + (gkyl_fend_fl - gkyl_fstart_fl) * \
            fstart_rans[i]
        fstart = np.round(fstart_fl)
        
        # Now create our Impurity object. 
        imp = PyImpurity(imp_atom_num, imp_mass, x, y, z, 
            imp_init_charge, fstart)
        
        if i > 0 and i % (num_imps / print_percent) == 0:
            print("  x = {:.4f}".format(x))
            print("  y = {:.4f}".format(y))
            print("  z = {:.4f}".format(z))
            print("  f = {:}".format(fstart))
        
        # Track transport of ions one frame at a time. Python indexes
        # the frames starting from 0, so we loop from 0 to gkyl_fend
        # minus fstart to account for that. 
        for f in range(0, gkyl_fend-fstart):
            
            # Add particle weight and count in our arrays. The particle
            # weights is just the timestep. You can think of this as
            # a particle in a cell is really only dt / (1 s) of a 
            # particle.
            imp_weight_arr[x_idx, y_idx, z_idx] += dt
            imp_count_arr[x_idx, y_idx, z_idx] += 1
            
            # Extract the plasma parameters at the current location.
            local_ne = gkyl_bkg["ne"][f, x_idx, y_idx, z_idx]
            local_te = gkyl_bkg["te"][f, x_idx, y_idx, z_idx]
            local_ni = gkyl_bkg["ni"][f, x_idx, y_idx, z_idx]
            local_ti = gkyl_bkg["ti"][f, x_idx, y_idx, z_idx]
            local_ex = gkyl_bkg["elecx"][f, x_idx, y_idx, z_idx]
            local_ey = gkyl_bkg["elecy"][f, x_idx, y_idx, z_idx]
            local_ez = gkyl_bkg["elecz"][f, x_idx, y_idx, z_idx]
            local_bz = gkyl_bkg["b"][x_idx, y_idx, z_idx]
            
            
            # Use ADAS to determine ionization/recombination 
            # probabilities, then pull a random number and see if we
            # are doing either of those. First load the rate 
            # coefficients.
            print("te={:.2f}  ti={:.2f}  ne={:.2e}  charge={}".format(
                local_te, local_ti, local_ne, imp.charge))
            iz_coef = oa.get_rate_coef(rate_df_iz, local_te, local_ne, 
                imp.charge)
            rc_coef = oa.get_rate_coef(rate_df_rc, local_te, local_ne, 
                imp.charge)
                
            # The probability of a ion ionizing/recombining during the 
            # timestep is then: 
            #   prob = coef * local_ne * dt
            iz_prob = iz_coef * local_ne * dt
            rc_prob = rc_coef * local_ne * dt
            
            # Pull random number, if it is less than the probability
            # then we have an event. Unfortunately we are pulling a
            # random number every loop, so this could probably be
            # optimized a bit.
            ran = np.random.random()
            if ran < iz_prob:
                imp.charge += 1
            #print("iz: ran={:.2e}  iz_prob={:.2e}".format(ran, iz_prob))
            ran = np.random.random()
            if ran < rc_prob:
                imp.charge -= 1
            #print("rc: ran={:.2e}  rc_prob={:.2e}".format(ran, rc_prob))
            
            
            # Check for wrong numbers?
            if imp.charge < 0:
                print("Error: imp_charge < 0 ({})".format(imp.charge))
            if imp.charge > imp.imp_atom_num:
                print("Error: imp_charge > {} ({})".format(imp.charge, 
                imp.imp_atom_num))
            if iz_prob > 1:
                print("Error: iz_prob ({}) > 1".format(iz_prob))
            if rc_prob > 1:
                print("Error: rc_prob ({}) > 1".format(rc_prob))
                
            # Perform step. imp object is updated within.
            print("Before imp_step: {:.6f}  {:.6f}  {:.6f}".format(
                imp.x, imp.y, imp.z))
            imp_step(imp, local_ex, local_ey, local_ez, local_bz, dt)
            print("After imp_step:  {:.6f}  {:.6f}  {:.6f}".format(
                imp.x, imp.y, imp.z))
            
            # Bounds checking. Treat x and z bounds as absorbing, and
            # the y bound as periodic.
            if imp.x <= gkyl_xmin or imp.x >= gkyl_xmax:
                break
            if imp.z <= gkyl_zmin or imp.z >= gkyl_zmax:
                break
            if imp.y <= gkyl_ymin:
                imp.y = gkyl_ymax + (imp.y - gkyl_ymin)
            elif imp.y >= gkyl_ymax:
                imp.y = gkyl_ymin + (imp.y - gkyl_ymax)
                
            
        
        # Free up memory.
        del imp

    # Calculate the impurity density. This is only for 3D cartesian
    # grid. Normalize to the sum then divide by the volume of each cell
    # to get values in m-3 (off by some unknown scalar).
    imp_dens_arr = imp_weight_arr / imp_weight_arr.sum()
    dx = gkyl_x[1] - gkyl_x[0]
    dy = gkyl_y[1] - gkyl_y[0]
    dz = gkyl_z[1] - gkyl_z[0]
    imp_dens_arr /= (dx * dy * dz)
    

    # Bundle up the things we care about and return.
    return_dict = {"imp_dens_arr": imp_dens_arr}
    return return_dict

def imp_step(imp, ex, ey, ez, bz, dt):
    """
    Update impurity location. We are using the GITR equations, defined
    in Younkin CPC 2021 (DOI: 10.1016/j.cpc.2021.107885):
    
    F = q(E + v x B) + Fc
    
    Fc = (see equation, too annoying to type out)
    """
    
    cdef float fx, fy, fz, dx, dy, dz, dvx, dvy, dvz, ev
    fx = 0.0
    fy = 0.0
    fz = 0.0
    dx = 0.0
    dy = 0.0
    dz = 0.0
    dvx = 0.0
    dvy = 0.0
    dvz = 0.0
    q = imp.charge * constants.ev
    
    # --------------------
    # Step in x direction.
    # --------------------
    # Electric field term
    fx += q * ex
    
    # Magnetic field term
    fx += q * imp.vy * bz
    
    # --------------------
    # Step in y direction.
    # --------------------
    # Electric field term
    fy += q * ey
    
    # Magnetic field term
    fy += -q * imp.vx * bz
    
    # --------------------
    # Step in z direction. 
    # --------------------
    # Electric field term
    fz += q * ez
    print ("ez = {:.2e}".format(ez))
    
    # Now calculate the step sizes.
    dvx = fx * dt / imp.mass
    dvy = fy * dt / imp.mass
    dvz = fz * dt / imp.mass
    #print("  dvx = {:.2e}  dvy = {:.2e}  dvz = {:.2e}".format(dvx, 
    #    dvy, dvz))
    dx = dvx * dt
    dy = dvy * dt
    dz = dvz * dt
    
    # Update imp object with new values.
    imp.x += dx
    imp.y += dy
    imp.z += dz
    imp.vx += dvx
    imp.vy += dvy
    imp.vz ++ dvz
    
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
