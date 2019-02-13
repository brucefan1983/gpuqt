# LSQT

An efficient implementation of linear scaling quantum transport (LSQT) methods which supports both pure CPU and GPU+CPU computations. This code can be used to obtain charge and spin transport properties of large systems described by a real-space tight-binding Hamiltonian. This is a work in progress. We aim to complete version 1.0 in a few months.

## References

The most original paper on this method is:
* [1] S. Roche and D. Mayou, Conductivity of Quasiperiodic Systems: A Numerical Study, Phys. Rev. Lett. 79, 2518 (1997). https://doi.org/10.1103/PhysRevLett.79.2518 

The major reference for the CUDA implementation is 
* [2] Z. Fan, A. Uppstu, T. Siro, and A. Harju, Efficient linear-scaling quantum transport calculations on graphics processing units and applications on electron transport in graphene, Comput. Phys. Commun. 185, 28 (2014). https://doi.org/10.1016/j.cpc.2013.08.009

This code was first published along with the following paper:
* [3] Z. Fan, V. Vierimaa, and Ari Harju, GPUQT: An efficient linear-scaling quantum transport code fully implemented on graphics processing units, Comput. Phys. Commun. 230, 113 (2018). https://doi.org/10.1016/j.cpc.2018.04.013

There is a comprehensive review article discussing the linear scaling quantum transport methods:
* [4] Zheyong Fan, Jose Hugo Garcia, Aron W. Cummings, Jose-Eduardo Barrios, Michel Panhans, Ari Harju, Frank Ortmann, and Stephan Roche, Linear Scaling Quantum Transport Methodologies, submitted to Reviews of Modern Physics, https://arxiv.org/abs/1811.07387

## File organizations

* After downloading and unpacking GPUQT, one can see two folders: `src` and `examples`. 

* The folder `src` contains all the source files (with suffix `.h` or `.cu`) of the main code and two makefiles. The source files are:
    * `main.cu`                            - the main function
    * `gpuqt.h` and `gpuqt.cu`             - the "driver function"
    * `sigma.h` and `sigma.cu`             - functions to obtain the transport properties
    * `model.h` and `model.cu`             - class to define the simulation model
    * `hamiltonian.h` and `hamiltonian.cu` - class to perform the matrix-related operations
    * `vector.h` and `vector.cu`           - class to perform the vector-related operations
    
* The folder `examples` contains two sub-folders with names `cpc2018` and `rmp`, containing examples in Refs. [3] and [4] mentioned above.

* There is also a file named `input.txt` in the `examples` folder, which is a "driver input file".

## Prerequisites

* The first version only supports GPU computing and it requires a computer/workstation/cluster equipped with one or more CUDA-enabled GPUs with compute capability of 2.0 or higher, a `g++` compiler and a `CUDA` toolkit. 
* Now the code also supports pure CPU computing and in this case it only requires a `g++` compiler.
* The code has only been tested in linux systems.

## Installing

* Go to `src` and 
    * type `make -f makefile.cpu` to build the CPU version. This will produce an executable called `lsqt_cpu` in the `src` folder.
    * type `make -f makefile.gpu` to build the GPU version. This will produce an executable called `lsqt_gpu` in the `src` folder.

## Running the examples

* Here I only show how to run the examples corresponding to Ref. [3].

* Go to the folders `examples/cpc2018/general/diffusive` and `examples/cpc2018/general/localized`, compile the `C++` codes with the `-std=c++11` option, and run the executables. This will produce all the input files (with suffix `.in`) for the examples.
  
* Go to the main folder where you can see the `src` folder and type `src/lsqt_gpu examples/input.txt` (or `src/lsqt_cpu examples/input.txt` if you use the CPU version) to run the examples. The data will be written into the output files (with suffix `.out`) in the `examples/cpc2018/general/diffusive` and `examples/cpc2018/general/localized` folders. If you run a simulation multiple times, new data will be appended to the existing output files.

## Analyzing the results

Go to the folders `examples/cpc2018/general/diffusive` and `examples/cpc2018/general/localized` and run the `MATLAB` scripts. You should get similar figures as in Ref. [3].


## Authors

* Zheyong Fan (Aalto University): Wrote the first working version of this code.

* Ville Vierimaa (Aalto University): Changed the code from the original C style to the current `C++` style and made many other improvements.

* Ari Harju (Aalto University): The supervisor of this project.

## Contact

* Zheyong Fan: brucenju(at)gmail.com; zheyong.fan(at)aalto.fi;

