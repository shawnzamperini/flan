# Flan
## Fully kinetic Monte Carlo to follow trace amount of impurities in a turbulent plasma. 

Flan is named after my dog, Flannery. It does not stand for anything. 

<img src="https://github.com/shawnzamperini/flan/blob/main/docs/flan_image.jpg" width="250">

Typical usage is to track a trace amount of impurities in a turbulent background generated by [Gkeyll](https://gkeyll.readthedocs.io/en/latest/install.html). Flan follows impurities in a plasma in a fully kinetic, Monte Carlo sense and is built on first-principles.

- The particles are tracked according to Lorentz force:
    $$\vec{F}=m(\vec{E} + \vec{v} \times \vec{B})$$
- Particle ionization and recombination is tracked via coupling to ADAS
- The collision model is based on the momentum slowing down time for a test particle as derived in the textbook *Plasma Physics and Fusion Energy* by Freidberg (see Ch. 9.3). Small changes are made to preserve the exact value of the Coulomb logarithm and to generalize it to particles colliding against any abritrary species. The relevant equations are:

    $$\nu_{12} = \frac{1}{4 \pi} \left( \frac{q_1^2 q_2^2 n_2 ln \Lambda}{\epsilon_0^2 m_1 m_r} \right) \left( \frac{1}{v_1^3 + 1.3 v_{T_2}^2} \right)$$

    Where generally 1 = impurity and 2 = the plasma species (electrons or deuterium ions). This frequency describes the rate at which an arbitrary impurity ion slows down, and can have a significant impact on the resulting impurity transport. For more details see the implementation in the [relevant source code](https://github.com/shawnzamperini/flan/blob/main/src/collisions.cpp).

A publication detailing Flan is anticipated in 2025. 

## Dependencies

CMake is used as a build system. The `mkdeps` shell script should handle installing all the dependencies, but they are: NetCDF, zlib, HDF5. Anaconda is also required, as the interface to Gkeyll (more specifically, the post-processing Gkeyll suite [postgkyl](https://github.com/ammarhakim/postgkyl/tree/main)) is written in python. A conda environment is automatically setup by `mkdeps` to handle this. 

## Installing Flan

Before installing flan, ensure you have a working installation of Anaconda on your machine. If you have that, then installing Flan should only require three commands. Navigate to the top directory of Flan and do the following:

1. Install dependencies: `machines/mkdeps.[machine].sh`
   - `machine` is where you are installing it on (probably `linux`)
   - This installs the dependencies in `$HOME/flansoft`. If you want them somewhere else, change the `FLANSOFT` variable in your `mkdeps` file
2. Generate build system: `cd build && cmake ..`
   - If you changed `FLANSOFT` above, modify the `cmake` command with `-DFLANSOFT=/path/to/flansoft`
3. Make Flan: `make`

This process installs the `flan` executable in the `build` directory and activates the `(flan)` conda environment. You can call `flan` using its full path or add an alias to its location, but it must be ran within the `(flan)` conda environment! 

## Regression Cases

To-do.

## Beginner's Guide

To-do.

## Post-processing

Interfacing with Flan data and making some basic plots are handled via a [python class interface](https://github.com/shawnzamperini/flan/blob/main/python/flan_plots.py). A tutorial showing how to use this will eventually be written, but most cases will only need `FlanPlots.plot_frame_xy` and/or `FlanPlots.plot_frames_xy`. Check the docstrings for more details on using those.  

## Upcoming Changes

- Investigate if GPUs can be used to speed up transport calculations
- Expand code to non-Cartesian geometries (e.g., tokamak coordinates)
