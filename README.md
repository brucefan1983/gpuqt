# gpuqt

A CUDA-implementation of a linear-scaling quantum transport method

## Installing

Go to "src" and type "make" to compile the gpuqt code. This will produce an executable "gpuqt" in the "src" folder.

## Running the examples

* Go to the folders "diffusive" and "localized", compile the C++ codes with the "-std=c++11" option, and run the executables.
* Go to the main folder where you can see the "src" folder and type "src/gpuqt examples/input.txt" to run the examples.

## Analyzing the results

Go to the folders "diffusive" and "localized" and run the matlab scripts. You should get similar figures as in the paper.
