========================
Installation and Testing 
========================

Dependencies
------------
CMake is used as a build system. The :literal:`mkdeps` shell script should handle installing all the dependencies, but they are: NetCDF, zlib, HDF5. Anaconda is also required, as the interface to Gkeyll (more specifically, the post-processing Gkeyll suite `postgkyl <https://github.com/ammarhakim/postgkyl/tree/main>`_) is written in python. A conda environment is automatically setup by :literal:`mkdeps` to handle this.


Installation
------------
Before installing Flan, ensure you have a working installation of Anaconda on your machine. If you are on a cluster, oftentimes you may need to load the conda module with something like :literal:`module load conda`, but ultimately each cluster has its own ways. If a system-wide installation is not available, you can install Anaconda in your home directory and run it from there. Once you have that, then installing Flan should only require three commands. Navigate to the top directory of Flan and do the following:

1. Install dependencies: :literal:`machines/mkdeps.[machine].sh` 

  - :literal:`machine` is where you are installing it on (probably linux)

  - This installs the dependencies in :literal:`$HOME/flansoft.` If you want them somewhere else, change the FLANSOFT variable in your mkdeps file

  - A conda environment named :literal:`flan` is created as part of this process, which you can activate with :literal:`conda activate flan`

2. Generate build system: :literal:`cd build && cmake ..`

  - If you changed :literal:`FLANSOFT` above, modify the cmake command with :literal:`-DFLANSOFT=/path/to/flansoft`

3. Make and install Flan: :literal:`make && make install`

This process installs the :literal:`flan` library in the :literal:`lib` directory and activates the :literal:`(flan)` conda environment. To run :literal:`flan`, one writes an input file in C++ and then compiles it with Flan. This process is detailed next, but it is important to note that the executable must be ran within the :literal:`(flan)` conda environment!


Testing
-------

CTest (part of CMake) is used to carry out tests. These are being written still. 


.. toctree::
