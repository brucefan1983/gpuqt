# GPUQT

A CUDA implementation of a linear-scaling quantum transport method. The major reference for the CUDA implementation is 

* Z. Fan, A. Uppstu, T. Siro, and A. Harju, Efficient linear-scaling quantum transport calculations on graphics processing units and applications on electron transport in graphene, Comput. Phys. Commun. 185, 28 (2014).

## File organizations

* After downloading and unpacking GPUQT, one can see two folders: "src" and "examples". 

* The folder "src" contains all the source files (with suffix .h or .cu) of the main code and a makefile. The source files are:
    * main.cu                          - the main function
    * common.h                         - some constants
    * gpuqt.h and gpuqt.cu             - the "driver function"
    * sigma.h and sigma.cu             - functions to obtain the transport properties
    * model.h and model.cu             - class to define the simulation model
    * hamiltonian.h and hamiltonian.cu - class to perform the matrix-related operations
    * vector.h and vector.cu           - class to perform the vector-related operations

* The folder "examples" contains two sub-folders with names "diffusive" and "localized", both containing the files "make_inputs.cpp" and "plot_results.m".

* There is also a file named "input.txt" in the "examples" folder, which is a "driver input file".

## Prerequisite

Need a CUDA-enabled GPU with compute capability of 2.0 or higher, a g++ compiler and a CUDA toolkit.

## Installing

Go to "src" and type "make" to compile the gpuqt code. This will produce an executable "gpuqt" in the "src" folder.

## Running the examples

* Go to the folders "diffusive" and "localized", compile the C++ codes with the "-std=c++11" option, and run the executables. This will produce all the input files (with suffix .in) for the examples.
  
* Go to the main folder where you can see the "src" folder and type "src/gpuqt examples/input.txt" to run the examples. The data will be written into the output files (with suffix .out) in the "diffusive" and "localized" folders. If you run a simulation multiple times, new data will be appended to the existing output files.

## Analyzing the results

Go to the folders "diffusive" and "localized" and run the MATLAB scripts. You should get similar figures as in the following paper:

* Z. Fan, V. Vierimaa, and Ari Harju, GPUQT: An efficient linear-scaling quantum transport code fully implemented on graphics processing units, arXiv:1705.01387 [physics.comp-ph], submitted to Comput. Phys. Commun. on May 3, 2017.

## Authors

* Zheyong Fan (Aalto University): wrote the first working version of this code.

* Ville Vierimaa (Aalto University): Changed the code from the original C stype to the current C++ style and made many other improvements.

* Ari Harju (Aalto University): The supervisor of this project.
