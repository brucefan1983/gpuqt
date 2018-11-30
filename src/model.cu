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

#include <fstream>
#include <sstream>
#include <chrono>




Model::Model(std::string input_dir)
{
    this->input_dir = input_dir;
    initialize_parameters();
    initialize_energy();
    if (requires_time)
        initialize_time();
    else
        time_step = 0;
    initialize_neighbor();
    initialize_positions();
    initialize_potential();
    initialize_hopping();
    random_state_real = new real[number_of_atoms]; 
    random_state_imag = new real[number_of_atoms];
 
    // Use higher accuracy clock for the RNG seed
    #ifdef DEBUG
        generator = std::mt19937(12345678);
    #else
       generator = std::mt19937
       (std::chrono::system_clock::now().time_since_epoch().count());
    #endif

    // We only need RNG for random phase generation 
    // so we may use interval [0, 2*PI] right away
    phase_distribution = std::uniform_real_distribution<real>(0, 2 * PI);
    std::cout << "Initialization complete\n" << std::endl;
}




Model::~Model()
{ 
    delete[] energy;
    delete[] time_step;
    delete[] potential;
    delete[] hopping_real;
    delete[] hopping_imag;
    delete[] neighbor_number;
    delete[] neighbor_list;
    delete[] xx;
    delete[] random_state_real;
    delete[] random_state_imag;
    delete[] x;
}




void Model::initialize_state(Vector& random_state)
{
    for (int n = 0; n < number_of_atoms; ++n)
    {  
        real random_phase = get_random_phase();
        random_state_real[n] = cos(random_phase);
        random_state_imag[n] = sin(random_phase);
    }
    random_state.copy_from_host(random_state_real, random_state_imag);
}




void Model::initialize_parameters()
{
    std::string filename = input_dir + "/para.in";
    std::cout << "\nReading " + filename << std::endl;
    std::ifstream input(filename);
    if (!input.is_open())
    {
        std::cout << "Error: cannot open " + filename << std::endl;
        exit(1);
    }
    std::string line;
    
    while (std::getline(input, line))
    {
        std::stringstream ss(line);
        std::string token;
        ss >> token;
        if (token == "") continue;
        if (token == "calculate_vac")
        {
            calculate_vac = true;
        }
        else if (token == "calculate_msd")
        {
            calculate_msd = true;
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
            std::cout << "Unknown identifier in " + input_dir + "/para.in:" << std::endl;
            std::cout << line << std::endl; 
        } 
    }
    input.close();
    
    if (calculate_vac || calculate_msd)
        requires_time = true;
    
    std::cout << "Finished reading " + filename << std::endl; 
    //Verify the used parameters
    std::cout << "- DOS will be calculated" << std::endl;
    if (calculate_vac)
        std::cout << "- VAC will be calculated" << std::endl;
    else
        std::cout << "- VAC is not calculated" << std::endl;
    if (calculate_msd)
        std::cout << "- MSD will be calculated" << std::endl;
    else
        std::cout << "- MSD is not calculated" << std::endl;    
    std::cout << "- Number of random vectors is " 
              << number_of_random_vectors << std::endl; 
    std::cout << "- Number of moments is " 
              << number_of_moments << std::endl;
    std::cout << "- Energy maximum is " << energy_max << std::endl; 
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
    input >> number_of_energy_points;
    energy = new real[number_of_energy_points];
    
    for (int n = 0; n < number_of_energy_points; ++n)
    {
        input >> energy[n];
    }
      
    input.close();
    std::cout << "Finished reading " + filename << std::endl;
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
    input >> number_of_steps_correlation;
    time_step = new real[number_of_steps_correlation];

    for (int n = 0; n < number_of_steps_correlation; ++n)
    {
        input >> time_step[n];
    }
    
    input.close();
    std::cout << "Finished reading " + filename << std::endl;         
}




void Model::initialize_neighbor()
{
    std::string filename = input_dir + "/neighbor.in";
    std::ifstream input(filename);
    
    if (!input.is_open()) 
    {
        std::cout <<"Error: cannot open " + filename << std::endl;
        exit(1);
    }
    input >> number_of_atoms >> max_neighbor;
    number_of_pairs = number_of_atoms * max_neighbor;
 
    neighbor_number = new int[number_of_atoms]; 
    neighbor_list = new int[number_of_pairs];

    for (int n = 0; n < number_of_atoms; ++n)
    {
        input >> neighbor_number[n];
        for (int m = 0; m < neighbor_number[n]; ++m)
        {        
            int index = n + m * number_of_atoms;
            input >> neighbor_list[index];
        }
    }
    input.close();
    std::cout << "Finished reading " + filename << std::endl; 
    std::cout << "- Number of atoms is " << number_of_atoms << std::endl;
    std::cout << "- Maximum neighbor number is " << max_neighbor << std::endl;            
}




real Model::get_random_phase()
{
    return phase_distribution(generator);
}




real reduce_distance(real d, real box)
{
    if (d > box/2.0)
        return d-box;
    if (d < -box/2.0)
        return d+box;
    else 
        return d;
}




void Model::initialize_positions()
{
    std::string filename = input_dir + "/position.in";
    std::ifstream input(filename);
    
    if (!input.is_open())
    {
        std::cout <<"Error: cannot open " + filename << std::endl;
        exit(1);
    }

    input >> box >> volume;   
    x = new real[number_of_atoms];

    for (int i=0; i<number_of_atoms; ++i)
        input >> x[i];
    input.close();
    std::cout << "Finished reading " + filename << std::endl;  
    std::cout << "- Box length along transport direction is " 
              << box << std::endl;
    std::cout << "- System volume is " << volume << std::endl;     
    std::cout << "- Calculating neighbor distances" << std::endl; 
  
    xx = new real[number_of_pairs];    
    for (int n = 0; n < number_of_atoms; ++n)
    {
        for (int m = 0; m < neighbor_number[n]; ++m)
        {        
            int index = n + m * number_of_atoms;
            xx[index] = reduce_distance(x[neighbor_list[index]] - x[n], box);
        }
    }
    std::cout << "- done" << std::endl;
}




void Model::initialize_potential()
{ 
    std::string filename = input_dir + "/potential.in";
    std::ifstream input(filename);
    bool nonzero_potential = true;
    if (!input.is_open())
    {
        std::cout <<"Could not open " + filename << std::endl;
        std::cout << "- Assuming zero onsite potential" << std::endl;
        nonzero_potential = false;
    }

    potential = new real[number_of_atoms];
    
    for (int n = 0; n < number_of_atoms; ++n)
    {
        if (nonzero_potential)
            input >> potential[n];
        else
            potential[n] = 0.0;
    }

    input.close();
    if (nonzero_potential)
        std::cout << "Finished reading " + filename << std::endl; 
}




void Model::initialize_hopping()
{
    std::string filename = input_dir + "/hopping.in";
    std::ifstream input(filename);

    /*
     type == 1 : complex hoppings
     type == 2 : real hoppings
     type == 3 : uniform hoppings (hoppings.in is not read)
    */
    int type = 0;
        
    if (!input.is_open())
    {
        std::cout <<"Could not open " + filename << std::endl;
        type = 3;
    }
    
    std::string first_line;
    

    if (type == 0)
        input >> first_line;
    else
        first_line = ".";
    
    if (first_line == "complex")
    {
        type = 1;
    }
    else if (first_line == "real")
    {
        type = 2;
    }
    else
    {
        type = 3;
        std::cout << "- Assuming uniform hoppings with strength 1" << std::endl;
    }
    
    hopping_real = new real[number_of_pairs]; 
    hopping_imag = new real[number_of_pairs];

    for (int n = 0; n < number_of_atoms; ++n)
    {
        for (int m = 0; m < neighbor_number[n]; ++m)
        {
            int index = n + m * number_of_atoms;
            if (type < 3)
                input >> hopping_real[index];
            else
                hopping_real[index] = 1.0;
            if (type == 1)
                input >> hopping_imag[index];
            else
                hopping_imag[index] = 0.0;
        }
    }
    input.close();
    if (type < 3)
        std::cout << "Finished reading " + filename << std::endl; 
    if (type == 1)
        std::cout << "- Hoppings had imaginary part" << std::endl;
    else if (type == 2)
        std::cout << "- Hoppings were real" << std::endl;
}




