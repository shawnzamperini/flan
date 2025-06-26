===================================
Welcome to the installation page!
===================================

Before installing flan, ensure you have a working installation of Anaconda on your machine. If you have that, then installing Flan should only require three commands. Navigate to the top directory of Flan and do the following:

Install dependencies: machines/mkdeps.[machine].sh
machine is where you are installing it on (probably linux)
This installs the dependencies in $HOME/flansoft. If you want them somewhere else, change the FLANSOFT variable in your mkdeps file
Generate build system: cd build && cmake ..
If you changed FLANSOFT above, modify the cmake command with -DFLANSOFT=/path/to/flansoft
Make and install Flan: make && make install
This process installs the flan library in the lib directory and activates the (flan) conda environment. To run flan, one writes an input file in C++ and then compiles it with Flan. This process is detailed next, but it is important to note that the executable must be ran within the (flan) conda environment!
