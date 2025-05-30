# Python interface to Gkeyll data via Postgkeyll. The high-level
# role of this script is to load the data using the well-maintained
# pgkyl package and then output in a format that is easily read
# in by Flan to do its calculations.
import postgkyl as pgkyl
import argparse
import sys
import numpy as np


def read_binary(args, comp, value_scale=1.0, bimax=False):
    """
    bimax (bool): Flag for if temperature and density are to be loaded from the
    BiMaxwellian files. The Gkeyll output is still changing, so in theory this
    would one day become unecessary, but until then we need to retain this
    capability to mesh with older sims.
    """
    
    # Go through one file at a time
    times = []
    grids = []
    values = []
    for frame in range(args.gkyl_frame_start, args.gkyl_frame_end + 1):

        # Assemble path and load
        if (args.gkyl_data_type == "density"):
            if bimax:
                data_dname = "BiMaxwellianMoments"
            else:
                data_dname = "M0"
        elif (args.gkyl_data_type == "temperature"):
            if bimax:
                data_dname = "BiMaxwellianMoments"
            else:
                data_dname = "prim_moms"
        elif (args.gkyl_data_type == "potential"):
            data_dname = "field"

        # bcart is actually the components of the unit vector of the magentic 
        # field in each X, Y, Z direction.
        elif (args.gkyl_data_type in ["magnetic_unit_X", "magnetic_unit_Y", 
            "magnetic_unit_Z"]):
            data_dname = "bcart"
        elif (args.gkyl_data_type == "jacobian"):
            data_dname = "jacobgeo"
        elif (args.gkyl_data_type == "magnetic_magnitude"):
            data_dname = "bmag"
        elif (args.gkyl_data_type in ["metric_coeff00", "metric_coeff01", 
            "metric_coeff02", "metric_coeff11", "metric_coeff12", 
            "metric_coeff22"]):
            data_dname = "gij"
        
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

        # Perform interpolation - add more options as I come across them
        if (data.ctx["basis_type"] == "serendipity"):
            interp_data = pgkyl.data.GInterpModal(data, data.ctx["poly_order"],
                "ms")
        else:
            print("Error! basis_type = {:} not supported yet."
                .format(data.ctx["basis_type"]))

        # Add time to time-dependent data.
        if (args.gkyl_data_type not in ["jacobian", "metric_coeff00", 
            "metric_coeff01", "metric_coeff02", "metric_coeff11", 
            "metric_coeff12", "metric_coeff22"]):
            times.append(data.ctx["time"])

        # Get grid and values. For each dimension, grid is one element larger
        # than values because it is the edges of each grid. Later in Flan we
        # can calculate cell centers, if needed. Note: We specify the component
        # here in interpolate.
        grid, value = interp_data.interpolate(comp)
        #grids.append(grid)

        # Option to multiply value by a scalar. This is useful when we need
        # to convert from velocity to temperature. 
        value *= value_scale
        values.append(value)

        # Ignore this comment below. We actually do want to do this so that the
        # magnetic field conforms to the other data arrays. One day if we
        # ever do electromagnetic simulations we won't need to change as much
        # code. 
        # If doing magnetic field, this is only to be done once since
        # we otherwise would be repeating the same information each frame.
        #if (args.gkyl_data_type == "magnetic_field"):
        #    break

    # Note: We intentionally just return grid and not grids since it is the
    # same for each time slice. If this changes, we can revisit but keep it 
    # simple for now.
    return times, grid, values

def save_csv(args, times, grid, values):
    """
    This functions creates three csv files for each array: times, grid, and 
    the values at each time.
    """
    
    # Use this as just a standard filename, no need to make it unique
    fname_base = "bkg_from_pgkyl_"

    # Add an informative header just for clarity's sake. The magnetic field
    # is just one array (electrostatic approximation for now), so there
    # are no times to write.
    if (args.gkyl_data_type not in ["magnetic_unit_X", "magnetic_unit_Y", 
        "magnetic_unit_Z", "magnetic_magnitude", "jacobian", "metric_coeff00", 
            "metric_coeff01", "metric_coeff02", "metric_coeff11", 
            "metric_coeff12", "metric_coeff22"]):
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
    with open(fname_base + args.gkyl_data_type + ".csv", "w") as f:
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
    Can optionally call this with standard arguments, they'll override anything
    passed in with the command line.
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

    # Parse arguments
    args = parser.parse_args()

    # Load binary file, specifying which components in the binary file has
    # the data we're after.
    if (args.gkyl_file_type == "binary"):

        # Load data that does not need to multiplied by anything and exists
        # in comp=0.
        if (args.gkyl_data_type in ["density", "potential", "jacobian", 
            "magnetic_magnitude"]):
            value_scale = 1.0
            comp = 0

        # Load temperature. This actually loads the velocity, which we convert
        # to temperature via T = value * m / q
        elif (args.gkyl_data_type == "temperature"):

            # Want to trigger a warning if the species mass was zero (probably
            # means it was accidentally not passed in). 
            if (args.gkyl_species_mass_amu == 0.0):
                print("Error! gkyl_species_mass_amu = 0.0. Did you forget" + \
                " to pass in this argument?")

            # Calculate mass of particle (kg) so we can convert to temperature
            m = 1.661e-27 * args.gkyl_species_mass_amu
            value_scale = m / 1.602e-19
            comp = 1

        # Magnetic field components
        elif (args.gkyl_data_type == "magnetic_unit_X"):
            value_scale = 1.0
            comp = 0
        elif (args.gkyl_data_type == "magnetic_unit_Y"):
            value_scale = 1.0
            comp = 1
        elif (args.gkyl_data_type == "magnetic_unit_Z"):
            value_scale = 1.0
            comp = 2

        # Metric coefficients. These make up a 3x3 array, so technically there 
        # are 9 different ones to load, but the matrix is symmetric so only
        # 6 are needed.
        elif (args.gkyl_data_type in ["metric_coeff00", "metric_coeff01", 
            "metric_coeff02", "metric_coeff11", "metric_coeff12", 
            "metric_coeff22"]):

            # No scaling needed
            value_scale = 1.0

            # Assign each specific component
            if args.gkyl_data_type == "metric_coeff00": comp = 0
            if args.gkyl_data_type == "metric_coeff01": comp = 1
            if args.gkyl_data_type == "metric_coeff02": comp = 2
            if args.gkyl_data_type == "metric_coeff11": comp = 3
            if args.gkyl_data_type == "metric_coeff12": comp = 4
            if args.gkyl_data_type == "metric_coeff22": comp = 5

        else:
            print("Error! Unrecognized gkyl_data_type: {}".format(args.gkyl_data_type))
            

        # When using BiMaxwellian files, the temperature needs an extra step.
        # this is because we have Tpar and Tperp. We want just the isotropic
        # temperature, which we can approximate as T = (Tpar + 2Tperp) / 3.
        # This is just the arithemetic mean. So we load both, and then estimate
        # it.
        if (args.gkyl_data_type == "temperature"):
            try:
                times, grid, values = read_binary(args, comp, value_scale)

            # Attempt looking for BiMaxwellian files then
            except FileNotFoundError:

                # Tpar/m = 2, Tperp/m = 3
                times, grid, values_par = read_binary(args, 2, value_scale, True)
                times, grid, values_perp = read_binary(args, 3, value_scale, True)
                values = (np.array(values_par) + 2 * np.array(values_perp)) / 3

        # Density could be stored in prim_moments or BiMaxwellian, try both.
        elif (args.gkyl_data_type == "density"):
            try: 
                times, grid, values = read_binary(args, comp, value_scale)
            except FileNotFoundError:

                # Is in comp = 0
                times, grid, values = read_binary(args, 0, value_scale, True)

        # Everything else has its own file
        else:
            times, grid, values = read_binary(args, comp, value_scale)

    else:
        print("Error! Only binary Gkyell files (.gkyl) are currently supported.")

    # Save to csv file with numpy.
    #print("Saving .csv files...")
    save_csv(args, times, grid, values)


if __name__ == "__main__":
    main()

