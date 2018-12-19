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




#pragma once
#include <random>
typedef double real;
class Vector;




class Model
{
public:
    Model(std::string input_dir);
    ~Model();
    void initialize_state(Vector& random_state);

    bool calculate_vac = false;
    bool calculate_msd = false;
    bool calculate_spin = false;

    int number_of_random_vectors = 1; 
    int number_of_atoms = 0; 
    int max_neighbor = 0;
    int number_of_pairs = 0;
    int number_of_energy_points = 0; 
    int number_of_moments = 1000; 
    int number_of_steps_correlation = 0;
    std::string input_dir;
    real energy_max = 10;

    real *energy;
    real *time_step;
    
    int *neighbor_number;
    int *neighbor_list;  
    real *xx;
    real *potential;
    real *hopping_real;
    real *hopping_imag;

    real volume;
    
private:

    void print_started_reading(std::string filename);
    void print_finished_reading(std::string filename);

    // for both lattice and general models
    void initialize_parameters();
    void initialize_energy();
    void initialize_time();

    // only for general model
    void initialize_neighbor();
    void initialize_positions();
    void initialize_potential();
    void initialize_hopping();
    void initialize_model_general();

    // only for lattice model
    void initialize_lattice_model();
    void add_anderson_disorder();
    void add_vacancies();
    void create_random_numbers(int, int, int*);
    void specify_vacancies(int*, int);
    void find_new_atom_index(int*, int*, int);
    
    bool requires_time = false;
    bool use_lattice_model = false;

    // disorder
    bool has_anderson_disorder = false;
    real anderson_disorder_strength;

    bool has_vacancy_disorder = false;
    int number_of_vacancies;
    
    std::mt19937 generator;
};




