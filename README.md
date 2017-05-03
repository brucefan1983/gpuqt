# GPUQT

A CUDA implementation of a linear-scaling quantum transport method. The major reference for the CUDA implementation is 

* Z. Fan, A. Uppstu, T. Siro, and A. Harju, Efficient linear-scaling quantum transport calculations on graphics processing units and applications on electron transport in graphene, Comput. Phys. Commun. 185, 28 (2014).

## Prerequisite

Need a CUDA-enabled GPU with compute capability of 2.0 or higher.

## Installing

Go to "src" and type "make" to compile the gpuqt code. This will produce an executable "gpuqt" in the "src" folder.

## Running the examples

* Go to the folders "diffusive" and "localized", compile the C++ codes with the "-std=c++11" option, and run the executables. This will produce all the input files (with suffix .in) for the examples.
  
* Go to the main folder where you can see the "src" folder and type "src/gpuqt examples/input.txt" to run the examples. The data will be written into the output files (with suffix .out) in the "diffusive" and "localized" folders. If you run a simulation multiple times, new data will be appended to the existing output files.

## Analyzing the results

Go to the folders "diffusive" and "localized" and run the MATLAB scripts. You should get similar figures as in the following paper (I will give the arXiv id when it's ready):

* Z. Fan, V. Vierimaa, and Ari Harju, GPUQT: An efficient linear-scaling quantum transport code fully implemented on graphics processing units

## Authors

* Zheyong Fan (Aalto University): wrote the first working version of this code.

* Ville Vierimaa (Aalto University): Changed the code from the original C stype to the current C++ style and made many other improvements.

* Ari Harju (Aalto University): The supervisor of this project.
