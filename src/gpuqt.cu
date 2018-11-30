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




#include "gpuqt.h"
#include "vector.h"
#include "hamiltonian.h"
#include "sigma.h"
#include "model.h"




void gpuqt(std::string input_directory)
{
    // Initialize model on the CPU
    Model model(input_directory);
    Hamiltonian H(model);
    Vector random_state(model);

    clock_t time_begin, time_finish;
    real time_used;
    
    // Loop over different random vectors
    for (int i = 0; i < model.number_of_random_vectors; ++i)
    {
        std::cout << "Starting vector number " << i+1 << std::endl;
        model.initialize_state(random_state);

        // Always calculate the DOS, since it is very cheap
        time_begin = clock(); 
        find_dos(model, H, random_state);
        time_finish = clock();
        time_used = real(time_finish - time_begin) / CLOCKS_PER_SEC;
        std::cout << "Time used for finding DOS = " 
                  << time_used << " s" << std::endl; 

        // Calculate the MSD only if you want to
        if (model.calculate_msd == 1)  
        {    
            time_begin = clock();
            find_msd(model, H, random_state);
            time_finish = clock();
            time_used = real(time_finish - time_begin) / CLOCKS_PER_SEC;
            std::cout << "Time used for finding MSD = " 
                      << time_used << " s" << std::endl;
        }

        // Calculate the VAC only if you want to
        if (model.calculate_vac == 1)  
        {
            time_begin = clock();
            find_vac(model, H, random_state);
            time_finish = clock();
            time_used = real(time_finish - time_begin) / CLOCKS_PER_SEC;
            std::cout << "Time used for finding VAC = " 
                      << time_used << " s" << std::endl;
        }
    }  
      	
}



