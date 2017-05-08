/*
    Copyright 2017 Zheyong Fan, Ville Vierimaa, and Ari Harju

    This file is part of GPUQT.

    GPUQT is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    GPUQT is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with GPUQT.  If not, see <http://www.gnu.org/licenses/>.
*/


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

