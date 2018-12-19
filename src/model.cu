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




#include "model.h"
#include "vector.h"
#include <random>
#include <iostream>
#include <fstream>
#include <sstream>
#include <chrono>

#define PI 3.141592653589793




Model::Model(std::string input_dir)
{
    // Use higher accuracy clock for the RNG seed
    #ifdef DEBUG
        generator = std::mt19937(12345678);
    #else
        generator = std::mt19937
        (std::chrono::system_clock::now().time_since_epoch().count());
    #endif

    // determine the input directory
    this->input_dir = input_dir;

    // read in para.in
    initialize_parameters();

    // read in energy.in and time_step.in
    initialize_energy();
    if (requires_time)
        initialize_time();
    else
        time_step = 0;

    // initialize the model
    if (use_lattice_model) // use a lattice model
    {
        initialize_lattice_model();
    }
    else // use general inputs to build the model
    {
        initialize_neighbor();
        initialize_positions();
        initialize_potential();
        initialize_hopping();
    }
}




Model::~Model()
{ 
    // other memory will be freed when constructing the Hamiltonian
    delete[] energy;
    delete[] time_step;
}




void Model::initialize_state(Vector& random_state)
{
    std::uniform_real_distribution<real> phase(0, 2 * PI);
    real *random_state_real = new real[number_of_atoms];
    real *random_state_imag = new real[number_of_atoms];

    // spin degeneracy is considered in perform_chebyshev_summation
    if (calculate_spin) // normalize to N/2 to remove spin degeneracy
    {
        for (int n = 0; n < number_of_atoms; n += 2)
        {
            real random_phase = phase(generator);
            random_state_real[n] = cos(random_phase);
            random_state_imag[n] = sin(random_phase);
            random_state_real[n+1] = 0.0;
            random_state_imag[n+1] = 0.0;
        }
    }
    else // normalize to N to keep spin degeneracy
    {
        for (int n = 0; n < number_of_atoms; ++n)
        {
            real random_phase = phase(generator);
            random_state_real[n] = cos(random_phase);
            random_state_imag[n] = sin(random_phase);
        }
    }
    random_state.copy_from_host(random_state_real, random_state_imag);
    delete[] random_state_real;
    delete[] random_state_imag;
}




void Model::print_started_reading(std::string filename)
{
    std::cout << std::endl;
    std::cout << "===========================================================";
    std::cout << std::endl;
    std::cout << "Started reading " + filename << std::endl;
    std::cout << std::endl;
}




void Model::print_finished_reading(std::string filename)
{
    std::cout << std::endl;
    std::cout << "Finished reading " + filename << std::endl;
    std::cout << "===========================================================";
    std::cout << std::endl << std::endl;
}




void Model::initialize_parameters()
{
    std::string filename = input_dir + "/para.in";
    std::ifstream input(filename);
    if (!input.is_open())
    {
        std::cout << "Error: cannot open " + filename << std::endl;
        exit(1);
    }
    print_started_reading(filename);

    std::string line;
    while (std::getline(input, line))
    {
        std::stringstream ss(line);
        std::string token;
        ss >> token;
        if (token == "") continue;
        if (token == "model")
        {
            ss >> use_lattice_model;
        }
        else if (token == "anderson_disorder")
        {
            has_anderson_disorder = true;
            ss >> anderson_disorder_strength;
        }
        else if (token == "vacancy_disorder")
        {
            has_vacancy_disorder = true;
            ss >> number_of_vacancies;
        }
        else if (token == "calculate_vac")
        {
            calculate_vac = true;
        }
        else if (token == "calculate_msd")
        {
            calculate_msd = true;
        }
        else if (token == "calculate_spin")
        {
            calculate_spin = true;
        }
        else if (token == "number_of_random_vectors")
        {
            ss >> number_of_random_vectors;
        }
        else if (token == "number_of_moments")
        {
            ss >> number_of_moments;
        }
        else if (token == "energy_max")
        {
            ss >> energy_max;
        }
        else
        {
            std::cout << "Unknown identifier in " + input_dir + "/para.in:" 
                      << std::endl;
            std::cout << line << std::endl;
        }
    }
    input.close();
    
    if (calculate_vac || calculate_msd)
        requires_time = true;
    
    //Verify the used parameters (make a seperate function later)
    if (use_lattice_model)
    {
        std::cout << "- Use lattice model" << std::endl;
        if (calculate_spin)
        {
            std::cout << "Error: lattice model does not support "
                      << "spin calculation yet" << std::endl;
            exit(1);
        }
    }
    else
        std::cout << "- Use general model" << std::endl;

    if (has_anderson_disorder)
    {
        std::cout << "- Add Anderson disorder with strength W = "
                  << anderson_disorder_strength << std::endl;
    }

    if (has_vacancy_disorder)
    {
        std::cout << "- Add " << number_of_vacancies
                  << " vacancies" << std::endl;
    }

    std::cout << "- DOS will be calculated" << std::endl;

    if (calculate_vac)
        std::cout << "- VAC will be calculated" << std::endl;
    else
        std::cout << "- VAC will not be calculated" << std::endl;

    if (calculate_msd)
        std::cout << "- MSD will be calculated" << std::endl;
    else
        std::cout << "- MSD will not be calculated" << std::endl;

    if (calculate_spin)
        std::cout << "- spin polarization will be calculated" << std::endl;
    else
        std::cout << "- spin polarization will not be calculated" << std::endl;

    if (calculate_spin && calculate_vac)
    {
        std::cout << "Error: spin and VAC cannot be calculated together"
                  << std::endl;
        exit(1);
    }

    if (calculate_spin && calculate_msd)
    {
        std::cout << "Error: spin and MSD cannot be calculated together"
                  << std::endl;
        exit(1);
    }

    std::cout << "- Number of random vectors is "
              << number_of_random_vectors << std::endl;
    std::cout << "- Number of moments is "
              << number_of_moments << std::endl;
    std::cout << "- Energy maximum is " << energy_max << std::endl;

    print_finished_reading(filename);
}




void Model::initialize_energy()
{
    std::string filename = input_dir + "/energy.in";
    std::ifstream input(filename);
    if (!input.is_open())
    {
        std::cout <<"Error: cannot open " + filename << std::endl;
        exit(1);
    }

    print_started_reading(filename);

    input >> number_of_energy_points;
    std::cout << "- number of energy points = "
              << number_of_energy_points 
              << std::endl;
    energy = new real[number_of_energy_points];

    for (int n = 0; n < number_of_energy_points; ++n)
    {
        input >> energy[n];
    }

    input.close();

    print_finished_reading(filename);
}




void Model::initialize_time()
{
    std::string filename = input_dir + "/time_step.in";
    std::ifstream input(filename);

    if (!input.is_open())
    {
        std::cout <<"Error: cannot open " + filename << std::endl;
        exit(1);
    }
    print_started_reading(filename);

    input >> number_of_steps_correlation;
    time_step = new real[number_of_steps_correlation];

    for (int n = 0; n < number_of_steps_correlation; ++n)
    {
        input >> time_step[n];
    }

    input.close();
    print_finished_reading(filename);
}




