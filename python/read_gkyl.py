# Python interface to Gkeyll data via Postgkeyll. The high-level
# role of this script is to load the data using the well-maintained
# pgkyl package and then output in a format that is easily read
# in by Flan to do its calculations.
import postgkyl as pgkyl
import argparse
import sys
import numpy as np


# Parameters that are time-independent
time_independent_data = ["jacobian", "metric_coeff00", 
    "metric_coeff01", "metric_coeff02", "metric_coeff11", "metric_coeff12", 
    "metric_coeff22"]
            
def load_binary_wrapper(args):
    """
    Wrapper for loading binary (.gkyl) data with postgkyl

    Parameters
    ----------
    args : Whatever parser.parse_args() return
        Contains all the arguments passed when calling read_gkyl.py

    Returns
    -------
    times : list
        List of all the times loaded

    grid : array
        3D (x,y,z) array of the grid coordinates

    values : array
        4D (t,x,y,z) array of the values loaded

    """
    time = None
    grid = None
    value = None

    # Get the name of file containing our data (data_dname), what component
    # in the file it's stored in (comp) and what it needs to be scaled by,
    # if at all (value_scale).
    data_dname, comp, value_scale = load_binary_params(args)

    # Go through one frame at a time, appending data as we go
    times = []
    values = []
    for frame in range(args.gkyl_frame_start, args.gkyl_frame_end + 1):

        # Assemble path
        path = load_binary_path(args, data_dname, frame)

        # If the user does not supply gkyl_moment_file_type, we can cycle
        # through the options until one doesn't fail.
        if (args.gkyl_moment_file_type is None and args.gkyl_data_type 
            in ["density", "temperature"]):
            for moment_file_type in ["m0m1m2", "maxwellian", "bimaxwellian"]:
                try:
                    
                    # Manually set moment_file_type here
                    args.gkyl_moment_file_type = moment_file_type

                    # Reload names and path and try again
                    data_dname, comp, value_scale = load_binary_params(args)
                    path = load_binary_path(args, data_dname, frame)
                    time, grid, value = load_binary(path, args, comp, value_scale)
                    break
                except FileNotFoundError:
                    pass

                #print("Error! Could not determine what type of file the " \
                #    "moments are stored in!")

        # If moment_file_type provided or not density/temperature
        else:
            time, grid, value = load_binary(path, args, comp, value_scale)

        # When loading temperature and the data is stored as BiMaxwellian, 
        # the data in value is actually T_par. An additional call is needed to 
        # get perpendicular temperature. Then T = (Tpar + 2Tperp) / 3.
        if (args.gkyl_data_type == "temperature" and 
            args.gkyl_moment_file_type == "bimaxwellian"):

            time, grid, value_perp = load_binary(path, args, 3, value_scale)
            value = (value + 2.0 * value_perp) / 3.0

        # Append. Don't need to do grid since it's the same every time. That
        # means grid just gets overwritten each loop, but whaddya gonna do
        # about it!
        times.append(time)
        values.append(value)

        # Some parameters don't have value for each frame, like the geometry
        # things, so we only need to do once and then we can break.
        if (args.gkyl_data_type in time_independent_data): break

    return times, grid, values

def load_binary_params(args):
    """
    Get the parameters needed to correctly load data from binary (.gkyl) file.

    Parameters
    ----------
    args : Whatever parser.parse_args() return
        Contains all the arguments passed when calling read_gkyl.py

    Returns
    -------
    data_dname : str
        The str used in constructing the path to the needed .gkyl file

    comp : int
        The component to load from the Gkeyll file (often 0)

    value_scale : float
        Multiplicative scale factor to convert loaded data to different
        units. Mainly used to convert the thermal velocity (vT^2) to eV for
        temperature.
    """

    data_dname = "null"
    comp = 0
    value_scale = 1.0

    # Density and temperature depend on how the moments were output from
    # Gkeyll. In all these, the thermal velocity is stored so one needs to
    # calculate the temperature (via the scale parameter). The options are:
    #   - M0M1M2: M0 = density, prim_moms = (?, vT^2, ?)
    #   - Maxwellian: (n, upar, vT^2)
    #   - BiMaxwellian: (n, upar, vT^2_par, vT^2_perp)
    # BiMaxwellian means we need to load twice, once for the parallel and
    # another for the perpendicular temperatures, 
    # then T = (Tpar + 2Tperp) / 3.
    if (args.gkyl_data_type == "temperature"):
        if (args.gkyl_moment_file_type == "m0m1m2"):
            data_dname = "prim_moms"
            comp = 1
        elif (args.gkyl_moment_file_type == "maxwellian"):
            data_dname = "MaxwellianMoments"
            comp = 2

        # Assigning to vT^2_par, will have additional load for vT^2_perp
        # (comp = 3) below to avoid a whole complicated logic since this is 
        # the only option that requires this.
        elif (args.gkyl_moment_file_type == "bimaxwellian"):
            data_dname = "BiMaxwellianMoments"
            comp = 2

        # Calculate mass of particle (kg) so we can convert to temperature
        m = 1.661e-27 * args.gkyl_species_mass_amu
        value_scale = m / 1.602e-19

    elif (args.gkyl_data_type == "density"):
        if (args.gkyl_moment_file_type == "m0m1m2"):
           data_dname = "M0"
        elif (args.gkyl_moment_file_type == "maxwellian"):
            data_dname = "MaxwellianMoments"
        elif (args.gkyl_moment_file_type == "bimaxwellian"):
            data_dname = "BiMaxwellianMoments"
        comp = 0
        value_scale = 1.0

    # Plasma potential
    elif (args.gkyl_data_type == "potential"):
        data_dname = "field"
        comp = 0
        value_scale = 1.0

    # bcart is the components of the unit vector of the magentic 
    # field in each X, Y, Z direction.
    elif (args.gkyl_data_type in ["magnetic_unit_X", "magnetic_unit_Y", 
        "magnetic_unit_Z"]):
        data_dname = "bcart"

        # Magnetic field components
        if (args.gkyl_data_type == "magnetic_unit_X"):
            comp = 0
        elif (args.gkyl_data_type == "magnetic_unit_Y"):
            comp = 1
        elif (args.gkyl_data_type == "magnetic_unit_Z"):
            comp = 2
        value_scale = 1.0

    # Jacobian
    elif (args.gkyl_data_type == "jacobian"):
        data_dname = "jacobgeo"
        comp = 0
        value_scale = 1.0

    # Magentic field strength
    elif (args.gkyl_data_type == "magnetic_magnitude"):
        data_dname = "bmag"
        comp = 0
        value_scale = 1.0

    # Metric coefficients
    elif (args.gkyl_data_type in ["metric_coeff00", "metric_coeff01", 
        "metric_coeff02", "metric_coeff11", "metric_coeff12", 
        "metric_coeff22"]):
        data_dname = "gij"

        # Assign correct component
        if args.gkyl_data_type == "metric_coeff00": comp = 0
        elif args.gkyl_data_type == "metric_coeff01": comp = 1
        elif args.gkyl_data_type == "metric_coeff02": comp = 2
        elif args.gkyl_data_type == "metric_coeff11": comp = 3
        elif args.gkyl_data_type == "metric_coeff12": comp = 4
        elif args.gkyl_data_type == "metric_coeff22": comp = 5

        value_scale = 1.0


    return data_dname, comp, value_scale

def load_binary_path(args, data_dname, frame):
    """
    Get full path to binary (.gkyl) file so it can be loaded

    Parameters
    ----------
    args : Whatever parser.parse_args() returns
        Contains all the arguments passed when calling read_gkyl.py

    data_dname : str
        The str used in constructing the path to the needed .gkyl file

    frame : int
        The frame number being loaded

    Returns
    -------
    path : str
        Full path to binary (.gkyl) file

    """

    # Some different ways each file name is constructed
    if (args.gkyl_data_type in ["density", "temperature"]):
        path = args.gkyl_dir + "/" + args.gkyl_case_name + "-" + \
            args.gkyl_species + "_" + data_dname + "_" + str(frame) \
            + ".gkyl"
    elif (args.gkyl_data_type == "potential"):
        path = args.gkyl_dir + "/" + args.gkyl_case_name + \
            "-" + data_dname + "_" + str(frame) \
            + ".gkyl"
    elif (args.gkyl_data_type in ["magnetic_unit_X", "magnetic_unit_Y", 
        "magnetic_unit_Z", "jacobian", "magnetic_magnitude", 
        "metric_coeff00", "metric_coeff01", "metric_coeff02", 
        "metric_coeff11", "metric_coeff12", "metric_coeff22"]):
        path = args.gkyl_dir + "/" + args.gkyl_case_name + \
            "-" + data_dname \
            + ".gkyl"

        # It seems the bmag file does not include this info, so it needs to
        # be passed in. This conditional statement allows this.
        if(args.gkyl_basis_type != "serendipity" or not  args.gkyl_poly_order):
            print("Error! Must pass in the interpolation basis type" + \
                " (gkyl_basis_type) and poly order (gkyl_poly_order)" + \
                " when saving the magnetic field data. Only " + \
                "serendipity is supported so far.")

    return path

def load_binary(path, args, comp, value_scale=1.0):
    """
    Load data from a binary (.gkyl) file

    Parameters
    ----------
    path : str
        The full path to the .gkyl file

    args : Whatever parser.parse_args() returns
        Contains all the arguments passed when calling read_gkyl.py

    comp : int
        The component to load from the Gkeyll file (often 0)

    value_scale : float
        Multiplicative scale factor to convert loaded data to different
        units. Mainly used to convert the thermal velocity (vT^2) to eV for
        temperature.

    Returns
    -------
    time : float
        The time loaded

    grid : array
        3D (x,y,z) array of the grid coordinates

    values : array
        3D (x,y,z) array of the values loaded

    """
    # Note: This accepts a comp argument but it is ignored for binary
    # filetypes. Very tricky. It is passed in below during interpolate. 
    data = pgkyl.data.GData(path)

    # Set the basis_type and poly_order if not included. It's really on
    # the user to ensure this matches among files until this information
    # is included in the bmag file.
    if (data.ctx["basis_type"] is None):
        data.ctx["basis_type"] = args.gkyl_basis_type
    if (data.ctx["poly_order"] is None):
        data.ctx["poly_order"] = args.gkyl_poly_order

    # Assign the values in data to the args. They should match. This allows
    # them to be written out later since we're using args to carry
    # everything around.
    args.gkyl_basis_type = data.ctx["basis_type"]
    args.gkyl_poly_order = data.ctx["poly_order"]

    # Perform interpolation - I think serendipity will be the only option
    if (data.ctx["basis_type"] == "serendipity"):
        interp_data = pgkyl.data.GInterpModal(data, data.ctx["poly_order"],
            "ms")
    else:
        print("Error! basis_type = {:} not supported yet."
            .format(data.ctx["basis_type"]))

    # Only time-dependent data has time in it
    if (args.gkyl_data_type not in ["jacobian", "metric_coeff00", 
        "metric_coeff01", "metric_coeff02", "metric_coeff11", 
        "metric_coeff12", "metric_coeff22"]):
        time = data.ctx["time"]
    else:
        time = []

    # Get grid and values. For each dimension, grid is one element larger
    # than values because it is the edges of each grid. Later in Flan we
    # can calculate cell centers, if needed. Note: We specify the component
    # here in interpolate.
    grid, value = interp_data.interpolate(comp)

    # Option to multiply value by a scalar. This is useful when we need
    # to convert from velocity to temperature. 
    value *= value_scale

    return time, grid, value

def save_csv(args, times, grid, values):
    """
    Creates three csv files for each array: times, grid, and the values at 
    each time. Parameters that do not depend on time skip the time step.

    Parameters
    ----------
    args : Whatever parser.parse_args() return
        Contains all the arguments passed when calling read_gkyl.py

    times : list
        List of all the times loaded

    grid : array
        3D (x,y,z) array of the grid coordinates

    values : array
        4D (t,x,y,z) array of the values loaded

    Returns
    -------
    N/A
    """
    
    # Use this as just a standard filename, no need to make it unique
    fname_base = "bkg_from_pgkyl_"

    # Add an informative header just for clarity's sake. The magnetic field
    # is just one array (electrostatic approximation for now), so there
    # are no times to write.
    if (args.gkyl_data_type not in time_independent_data):
        header = (
                 "# The times of each Gkeyll frame. The first line is an\n"
                 "# integer telling how many values follow.\n"
                 )
                 
        with open(fname_base + "times.csv", "w") as f:
            f.write(header)
            num_vals = "{:d}\n".format(len(times))
            f.write(num_vals)
            np.savetxt(f, times)

    # Add an informative header just for clarity's sake
    header = (
             "# The grid loaded by postgkyl. These are the nodes of the grid,\n"
             "# so there is one additional value in each direction relative to\n" 
             "# the values array. There are three arrays here, one for each\n"
             "# direction (x, y, z). The format of this file is a line with\n"
             "# three integers telling how many values each array is, and then\n"
             "# all the arrays.\n"
             )             

    # Save the grid file.
    with open(fname_base + "grid.csv", "w") as f:
        f.write(header)
        num_vals = "{:d} {:d} {:d}\n".format(len(grid[0]), len(grid[1]), len(grid[2]))
        f.write(num_vals)
        for i in range(0, 3):
            np.savetxt(f, grid[i])

    # Add an informative header just for clarity's sake
    header = ( 
             "# The values for the specified type of data requested from\n"
             "# Gkeyll. The first row are integers: the number of frames,\n"
             "# and then the shape of each array (x, y, z dimensions).\n"
             "# Then each array for each frame is printed out flattened\n"
             "# using C-style indexing.\n"
             ) 
    
    # Data like density and temperature could be ion or electron (or some other
    # species), so we need to add that extra qualifier.
    if args.gkyl_data_type in ["density", "temperature"]:
        fname = fname_base + args.gkyl_species + "_" + \
            args.gkyl_data_type + ".csv"
    else:
        fname = fname_base + args.gkyl_data_type + ".csv"

    with open(fname, "w") as f:
        f.write(header)
        num_vals = "{:d} {:d} {:d} {:d}\n".format(len(values), values[0].shape[0], 
                values[0].shape[1], values[0].shape[2]) 
        f.write(num_vals)
        for i in range(0, len(values)):
            np.savetxt(f, values[i].flatten())

    # Save the interpolation settings
    header = (
            "# The DG interpolation options used. These should probably\n"
            "# be the same for each data type in a run, so in theory this\n"
            "# will be overwriting a file with the same data each time.\n"
            "# The first line is the basis type, second is the poly order.\n"
            )
    with open(fname_base + "interp_settings.csv", "w") as f:
        f.write(header)
        f.write(args.gkyl_basis_type + "\n")
        f.write(str(args.gkyl_poly_order) + "\n")


def main():
    """
    Command line utility to load data from Gkeyll using postgkyl and save the
    data in .csv files to be read in by Flan. Effectively the coupling between
    Gkeyll output and Flan input is here.
    """

    # Valid options for gkyl_data_type
    gkyl_data_type_opts = ["density", "temperature", "potential", "magnetic_unit_X", 
        "magnetic_unit_Y", "magnetic_unit_Z", "jacobian", "magnetic_magnitude", 
        "metric_coeff00", "metric_coeff01", "metric_coeff02", "metric_coeff11",
        "metric_coeff12", "metric_coeff22"]

    # Parse command line arguments
    parser = argparse.ArgumentParser(
        description="Read Gkeyll data and save it for Flan to read in")

    # Command line arguments
    parser.add_argument("--gkyl_dir", type=str, 
        help="full path to the directory containing Gkyell simulation results")
    parser.add_argument("--gkyl_case_name", type=str, 
        help="name of Gkyell simulation to load")
    parser.add_argument("--gkyl_file_type", type=str, 
        help="type of file for the Gkyell data", 
        choices=["adios1", "adios2", "binary"])
    parser.add_argument("--gkyl_frame_start", type=int, 
        help="first Gkyell frame to read in")
    parser.add_argument("--gkyl_frame_end", type=int, 
        help="last Gkyell frame to read in")
    parser.add_argument("--gkyl_species", type=str, 
        help="name of the species to load")
    parser.add_argument("--gkyl_species_mass_amu", type=float, 
        help="mass of the species in amu")
    parser.add_argument("--gkyl_data_type", type=str, 
        help="name of the type of data file to load", 
        choices=gkyl_data_type_opts)
    parser.add_argument("--gkyl_basis_type", type=str, help="DG basis type", 
        choices=["serendipity"])
    parser.add_argument("--gkyl_poly_order", type=int, help="DG poly order")
    parser.add_argument("--gkyl_moment_file_type", type=str, 
        help="type of moment files containing density and temperature",
        choices=["m0m1m2", "maxwellian", "bimaxwellian"])

    # Parse arguments
    args = parser.parse_args()

    # Load binary file, specifying which components in the binary file has
    # the data we're after.
    if (args.gkyl_file_type == "binary"):
        times, grid, values = load_binary_wrapper(args)
    else:
        print("Error! Only binary Gkyell files (.gkyl) are currently supported.")

    # Save to csv file with numpy.
    save_csv(args, times, grid, values)


if __name__ == "__main__":
    main()

