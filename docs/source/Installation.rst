===================================
Welcome to the dependencies/installation page!
===================================
Dependencies
------------
CMake is used as a build system. The :literal:`mkdeps` shell script should handle installing all the dependencies, but they are: NetCDF, zlib, HDF5. Anaconda is also required, as the interface to Gkeyll (more specifically, the post-processing Gkeyll suite `postgkyl <https://github.com/ammarhakim/postgkyl/tree/main>`_) is written in python. A conda environment is automatically setup by :literal:`mkdeps` to handle this.


Installation
-------------
Before installing flan, ensure you have a working installation of Anaconda on your machine. If you have that, then installing Flan should only require three commands. Navigate to the top directory of Flan and do the following:

1. Install dependencies: :literal:`machines/mkdeps.[machine].sh` 

    -:literal:`machine` is where you are installing it on (probably linux)

    -This installs the dependencies in :literal:`$HOME/flansoft.` If you want them somewhere else, change the FLANSOFT variable in your mkdeps file

2. Generate build system: :literal:`cd build && cmake ..`

  -If you changed :literal:`FLANSOFT` above, modify the cmake command with :literal:`-DFLANSOFT=/path/to/flansoft`

3. Make and install Flan: :literal:`make && make install`

This process installs the :literal:`flan` library in the :literal:`lib` directory and activates the :literal:`(flan)` conda environment. To run :literal:`flan`, one writes an input file in C++ and then compiles it with Flan. This process is detailed next, but it is important to note that the executable must be ran within the :literal:`(flan)` conda environment!


