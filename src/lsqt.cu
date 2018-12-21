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




#include "lsqt.h"
#include "vector.h"
#include "hamiltonian.h"
#include "sigma.h"
#include "model.h"
#include <iostream>
typedef double real;




static void print_started_random_vector(int i)
{
    std::cout << std::endl;
    std::cout << "===========================================================";
    std::cout << std::endl;
    std::cout << "Started  simulation with random vector number " 
              << i << std::endl;
    std::cout << std::endl;
}



static void print_finished_random_vector(int i)
{
    std::cout << std::endl;
    std::cout << "Finished simulation with random vector number " 
              << i << std::endl; 
    std::cout << "===========================================================";
    std::cout << std::endl << std::endl;
}




void lsqt(std::string input_directory)
{
    // Initialize model on the CPU
    Model model(input_directory);
    Hamiltonian H(model);
    Vector random_state(model.number_of_atoms);

    clock_t time_begin, time_finish;
    real time_used;

    // Loop over different random vectors
    for (int i = 0; i < model.number_of_random_vectors; ++i)
    {
        print_started_random_vector(i);

        int orbital = -1; // using random vectors rather than a local orbital
        model.initialize_state(random_state, orbital);

        // Always calculate the DOS, since it is very cheap
        time_begin = clock(); 
        find_dos(model, H, random_state, 0);
        time_finish = clock();
        time_used = real(time_finish - time_begin) / CLOCKS_PER_SEC;
        std::cout << "- Time used for finding DOS = " 
                  << time_used << " s" << std::endl; 

        // Calculate the MSD only if you want to
        if (model.calculate_msd == 1)  
        {    
            time_begin = clock();
            find_msd(model, H, random_state);
            time_finish = clock();
            time_used = real(time_finish - time_begin) / CLOCKS_PER_SEC;
            std::cout << "- Time used for finding MSD = " 
                      << time_used << " s" << std::endl;
        }

        // Calculate the VAC only if you want to
        if (model.calculate_vac == 1)  
        {
            time_begin = clock();
            find_vac(model, H, random_state);
            time_finish = clock();
            time_used = real(time_finish - time_begin) / CLOCKS_PER_SEC;
            std::cout << "- Time used for finding VAC = " 
                      << time_used << " s" << std::endl;
        }

        // Calculate the spin polarization only if you want to
        if (model.calculate_spin == 1)  
        {
            time_begin = clock();
            find_spin_polarization(model, H, random_state);
            time_finish = clock();
            time_used = real(time_finish - time_begin) / CLOCKS_PER_SEC;
            std::cout << "- Time used for finding spin polarization = " 
                      << time_used << " s" << std::endl;
        }

        print_finished_random_vector(i);
    }

    // Calculate the LDOS only if you want to
    if (model.calculate_ldos)
    {
        time_begin = clock();
        // loop over the local orbitals
        for (int i = 0; i < model.number_of_local_orbitals; ++i)
        {
            int orbital = model.local_orbitals[i];
            model.initialize_state(random_state, orbital);
            find_dos(model, H, random_state, orbital);
            std::cout << "- Finished orbital " << orbital << std::endl;
        }
        time_finish = clock();
        time_used = real(time_finish - time_begin) / CLOCKS_PER_SEC;
        std::cout << "- Time used for finding LDOS = "
                  << time_used << " s" << std::endl;
    }
}



