# distutils: language=c++

# Copying the example found here:
# https://cython.readthedocs.io/en/stable/src/userguide/wrapping_CPlusPlus.html#wrapping-cplusplus
# Main particle following loop. Must be compiled with as cython code. 
# See setup.py file for details. 
import sys
import numpy as np
import constants
import time
from tqdm import tqdm
from scipy.interpolate import RegularGridInterpolator
import cython

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

def follow_impurities(input_opts, gkyl_bkg):
    """
    This function contains the primary impurity ion following routine.
    """
    
    # Define variables to use as much C++ as possible here.
    cdef int i, j, imp_zstart_int, x_idx, y_idx, z_idx, fstart, fend
    cdef int loop_counts, iz_warn_count, rc_warn_count
    cdef float x, y, z, local_ne, local_te, local_ti, local_ni
    cdef float local_ex, local_ey, local_ez, iz_coef, rc_coef, iz_prob
    cdef float rc_prob, ran, local_bz, print_percent, perc_done
    cdef float avg_time_followed, local_dtedx, local_dtedy, local_dtedz
    cdef float local_dtidx, local_dtidy, local_dtidz, local_viz
    cdef float local_stop_time, fstart_fl, fend_fl
    
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
    cdef bint imp_xyz_tracks = input_opts["imp_xyz_tracks"]
    
    # If we interpolated frames the end frame for our context will be
    # larger by gkyl_num_interp + 1. See read_gkyl.interp_frames for
    # more info. 
    #gkyl_num_interp = input_opts["gkyl_num_interp"]
    #if gkyl_num_interp is not None:
    #    gkyl_fend *= (gkyl_num_interp + 1)
    
    
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
        
    elif imp_zstart_opt == "single_value":
        imp_zstart_int = 2
        
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
    
    # Load the Gkeyll background into memory views. 
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
    
    # Calculate the stopping time (Stangeby Eq. 6.35) ahead of time.
    # Calculate as python numpy array then assign to C array.
    ln_alpha = 15.0
    ion_amu = 2.04  # Hardcoding deuterium in.
    tmp_stop_time = 1.47e13 * imp_amu * gkyl_bkg["ti"] * \
        np.sqrt(gkyl_bkg["ti"] / ion_amu) / \
        ((1 + ion_amu / imp_amu) * gkyl_bkg["ni"] *  
        np.power(imp_atom_num, 2) * ln_alpha)
    cdef float[:,:,:,:] stop_time = tmp_stop_time
    
    # Calculate the alpha and beta parameters that go in front of the 
    # electron and ion temperature gradient forces, respectively. 
    mu = imp_mass / (imp_mass + constants.mi)
    cdef float alpha = 0.71 * np.power(imp_atom_num, 2)
    cdef float beta = 3 * (mu + 5 * np.sqrt(2) * 
        np.square(imp_atom_num) * (1.1 * np.power(mu, 5/2) - 0.35 * 
        np.power(mu, 3/2)) - 1) / (2.6 - 2 * mu + 5.4 * np.square(mu))
    
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
    
    # Remember, these lists are 0-index, so for a specific charge number
    # we will want to index it with charge-1. The functions are called
    # with interps[charge-1](local_te, local_ne).
    interps_rc = create_interps(rate_df_rc, imp_atom_num)
    interps_iz = create_interps(rate_df_iz, imp_atom_num)
    
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
    
    # Lists to record the impurity ion tracks. I am going to hardcode
    # a safegaurd in here for now, since I don't want to accidentaly run
    # this and make a massive output file. I think you shoul dbe able to
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
    
    # Initialize some loop variables.
    loop_counts = 0
    iz_warn_count = 0
    rc_warn_count = 0
    avg_time_followed = 0.0
    imp_foll_tstart = time.time()

    # Loop through tracking the full path of one impurity at a time.
    print("Beginning impurity following...")
    #for i in tqdm(range(0, num_imps)):
    for i in range(num_imps):
        
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
            z_idx = np.argmin(np.abs(zstart_cdf[x_idx, y_idx] 
            - zstart_rans[i]))
            z = gkyl_z[z_idx]
        elif imp_zstart_int == 2:
            z = imp_zstart_val
        
        # Assume impurity can start at any frame (this means we are 
        # assuming a constant impurity source). This could lead to 
        # low statistics at the starting frames, should be studied. 
        #fstart_fl = gkyl_fstart_fl + (gkyl_fend_fl - gkyl_fstart_fl) * \
        #    fstart_rans[i]
        #fstart = np.round(fstart_fl)
        fstart_fl = fend_fl * fstart_rans[i]
        fstart = np.round(fstart_fl)
        
        # Now create our Impurity object. 
        imp = PyImpurity(imp_atom_num, imp_mass, x, y, z, 
            imp_init_charge, fstart)
        
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
        #for f in range(fstart-gkyl_fstart, gkyl_fend-gkyl_fstart):
        for f in range(fstart, fend):
            
            # Update nearest index values.
            x_idx = np.argmin(np.abs(imp.x - gkyl_x))
            y_idx = np.argmin(np.abs(imp.y - gkyl_y))
            z_idx = np.argmin(np.abs(imp.z - gkyl_z))
            
            # Add particle weight and count in our arrays. The particle
            # weights is just the timestep. You can think of this as
            # a particle in a cell is really only dt / (1 s) of a 
            # particle.
            imp_weight_arr[f, x_idx, y_idx, z_idx] += dt
            imp_count_arr[f, x_idx, y_idx, z_idx] += 1
            
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
            else:
                local_dtedx = 0.0
                local_dtedy = 0.0
                local_dtedz = 0.0
                local_dtidx = 0.0
                local_dtidy = 0.0
                local_dtidz = 0.0
            
            # Use ADAS to determine ionization/recombination 
            # probabilities, then pull a random number and see if we
            # are doing either of those. First load the rate 
            # coefficients by calling the function for this charge state
            # that we already loaded before the main loop.
            #print("te={:.2f}  ti={:.2f}  ne={:.2e}  charge={}".format(
            #    local_te, local_ti, local_ne, imp.charge))
            rc_coef = interps_rc[imp.charge-1]((local_te, local_ne))
            iz_coef = interps_iz[imp.charge-1]((local_te, local_ne))
                
            # The probability of a ion ionizing/recombining during the 
            # timestep is then: 
            #   prob = coef [m3/s] * local_ne [m-3] * dt [s]
            iz_prob = iz_coef * local_ne * dt
            rc_prob = rc_coef * local_ne * dt
            
            # Pull random number, if it is less than the probability
            # then we have an event. Unfortunately we are pulling a
            # random number every loop, so this could probably be
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
                print("Error: imp_charge < 0 ({})".format(imp.charge))
            if imp.charge > imp.imp_atom_num:
                print("Error: imp_charge > {} ({})".format(imp.charge, 
                imp.imp_atom_num))
            if iz_prob > 1:
                #print("Error: iz_prob ({}) > 1".format(iz_prob))
                iz_warn_count += 1
            if rc_prob > 1:
                #print("Error: rc_prob ({}) > 1".format(rc_prob))
                rc_warn_count += 1
                
            # Perform step. imp object is updated within.
            #print("Before imp_step: {:.6f}  {:.6f}  {:.6f} ({:.3e} {:.3e}  {:.3e})".format(imp.x, imp.y, imp.z, imp.vx, 
            #    imp.vy, imp.vz))
            imp_step(imp, local_ex, local_ey, local_ez, local_bz, dt, 
                local_dtedx, local_dtedy, local_dtedz, local_dtidx, 
                local_dtidy, local_dtidz, local_viz, local_stop_time, 
                alpha, beta, coll_forces)
            #print("After imp_step: {:.6f}  {:.6f}  {:.6f} ({:.3e} {:.3e}  {:.3e})".format(imp.x, imp.y, imp.z, imp.vx, 
            #    imp.vy, imp.vz))
            #print("{}-{}: x={:.4f}  y={:.4f}  z={:.4f}  ({}, {}, {})"
            #    .format(i, f, imp.x, imp.y, imp.z, x_idx, y_idx, z_idx))
            
            loop_counts += 1
            avg_time_followed += dt
            
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
            
        # Free up memory (if necessary since it gets overwritten later).
        del imp

    # Calculate how long the main loop took.
    imp_foll_tend = time.time()
    cpu_time_used = imp_foll_tend - imp_foll_tstart
    print("Time spent following impurities: {:.1f} s".format(cpu_time_used))

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
    
    # Average time following an impurity ion.
    avg_time_followed /= num_imps
    
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

    # Bundle up the things we care about and return.
    return_dict = {
        "imp_dens_arr": imp_dens_arr, 
        "dt": dt, 
        "avg_time_followed": avg_time_followed, 
        "cpu_time_used": cpu_time_used
        }
        
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

