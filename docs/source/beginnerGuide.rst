=========================================================================================
Beginner's Guide
=========================================================================================
It is useful to create a directory to contain all your :literal:`flan` cases, it can be anywhere (don't put it in the repository directory, bad practice). For this section we will assume it is called :literal:`flandir`. Within :literal:`flandir`, you can use the following script to setup a simulation directory named :literal:`test`:

.. code-block:: bash

   (flan) $ /path/to/flan/scripts/flan_setup.sh test

This will create a directory called :literal:`test` and place the needed files in it. The input file is :literal:`test.cpp`. If you open this up, you will notice a function :literal:`mapc2p`, a function called :literal:`create_inputs` and :literal:`main`. For now, all you need to look at is :literal:`create_inputs`. This is where all the input options are entered. You will notice some input options are provided to get you started. It is good to be familiar with the Gkeyll simulation you are running :literal:`flan` with so that you can tell it where the impurities start. A possible bare-minimum set of input options could look like this:

.. code-block:: bash

  inpts["case_name"] = "test";
  inpts["imp_num"] = 10000;
  
  // Tell Flan where to find the Gkeyll files and how many to use. This is one
  // that used the simple helical approximation, so the coordinates are already
  // Cartesian (this means nothing needs to be done with mapc2p).
  inpts["gkyl_dir"] = "/home/zamp/gkyldir/d3d-167196-v6-gpu";
  inpts["gkyl_casename"] = "d3d-167196-v6-gpu";
  inpts["gkyl_frame_start"] = 600;
  inpts["gkyl_frame_end"] = 699;
  inpts["gkyl_elec_name"] = "elc";
  inpts["gkyl_ion_name"] = "ion";
  
  // Setup simulation to follow W ions
  inpts["imp_num"] = 10000
  inpts["imp_mass_amu"] = 183.84;
  
  // I know from this Gkeyll simulation that this corresponds to the "left"
  // or inner edge of the simulation, so start W ions there and we will
  // watch them transport from there.
  inpts["imp_xmin"] = 2.315;
  inpts["imp_xmax"] = 2.315;
  
  // Randomly start the ions anywhere in the y direction
  inpts["imp_ystart_opt"] = "range";
  
  // Need to know where to find the ADAS files and what year to use
  inpts["openadas_root"] = "/home/zamp/flandir/openadas";
  inpts["openadas_year"] = 50;

You can leave the rest of the input file alone for now. To run the simulation, you must compile it first to create an executable, then run the executable. This is easily done (don't forget to make sure the :literal:`flan` conda environment is active):

.. code-block:: bash

  (flan) $ make
  (flan) $ ./test

Once finished, we can plot the data using the provided :literal:`flan_plots` python plotting library. Note: You must have the :literal:`flan` conda environment active

.. code-block:: bash

  (flan) $ ipython
  In [1]: import flan_plots
  In [2]: fp = flan_plots.FlanPlots("test.nc")
  In [3]: fp.plot_frames_xy("imp_density", 0, 99, 0.0, animate_cbar=True, vmin=1e-6, vmax=1e-3, norm_type="log", xlabel="R-Rsep (m)", ylabel="Binormal (m)", cbar_label="W Density (arb.)")

This creates an animated plot, a screenshot of which is shown below.

.. image:: flan/docs/flan_ex_v1.png

Data can be extracted for more detailed analysis. This will be covered in future :literal:`flan_plots` documentation (one day).

