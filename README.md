# GPUQT

A CUDA-implementation of a linear-scaling quantum transport method. The major reference for the GPU implementation is 

* Z. Fan, A. Uppstu, T. Siro, and A. Harju, Efficient linear-scaling quantum transport calculations on graphics processing units and applications on electron transport in graphene, Comput. Phys. Commun. 185, 28 (2014).

## Prerequisite

Need a CUDA-enabled GPU with compute capability of 2.0 or higher.

## Installing

Go to "src" and type "make" to compile the gpuqt code. This will produce an executable "gpuqt" in the "src" folder.

## Running the examples

* Go to the folders "diffusive" and "localized", compile the C++ codes with the "-std=c++11" option, and run the executables.
  This will produce all the input files for the examples.
  
* Go to the main folder where you can see the "src" folder and type "src/gpuqt examples/input.txt" to run the examples.

## Analyzing the results

Go to the folders "diffusive" and "localized" and run the matlab scripts. You should get similar figures as in the paper (I will give the arXiv id when it's ready).
