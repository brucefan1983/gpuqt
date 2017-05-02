/*----------------------------------------------------------------------------80
    This is the main function of GPUQT (GPU Quantum Transport).

    Authors: Zheyong Fan <brucenju@gmail.com> <zheyongfan@163.com> 
             Ville Vierimaa <ville.v.vierimaa@aalto.fi>
------------------------------------------------------------------------------*/

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <vector>
#include <iostream>
#include <fstream>
#include <string>

#include "common.h"
#include "gpuqt.h"

//----------------------------------------------------------------------------80
int main(int argc, char *argv[])
{	
    if (argc != 2)
    {
        std::cout << "Usage: src/gpuqt input.txt" << std::endl;
        exit(1);
    }
	
    std::ifstream input(argv[1]);
    if (!input.is_open())
    {
        std::cout << "Failed to open " << argv[1] << std::endl;
        exit(1);
    }		

    std::string directory;
    while (std::getline(input, directory))
    {
 		if (directory == "")
 			continue;
        std::cout << std::endl;
        std::cout << "===========================================" << std::endl;
        std::cout << "Run KGQT simulation for " << directory << std::endl; 
        std::cout << "===========================================" << std::endl;

        clock_t time_begin = clock();
         
        // call the driver function
        gpuqt(directory);

        clock_t time_finish = clock();
        real time_used = real(time_finish - time_begin) / CLOCKS_PER_SEC;

        std::cout << std::endl;
        std::cout << "===========================================" << std::endl;
        std::cout << "Total time used for " << directory << " = " << time_used <<" s" << std::endl; 
        std::cout << "===========================================" << std::endl;
    }

    return 0;
}

