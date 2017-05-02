#pragma once
#include "common.h"
#include <random>

class Vector;

class Model
{
public:
	Model(std::string input_dir);
	~Model();
	void initialize_state(Vector& random_state);

    bool calculate_vac = false;
    bool calculate_msd = false;

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
    void initialize_parameters();
    void initialize_energy();
    void initialize_time();
    void initialize_neighbor();
    void initialize_positions();
    void initialize_potential();
    void initialize_hopping();
	
	real get_random_phase();

    real *random_state_real;
    real *random_state_imag;    	
    
    real* x;
    real box;
    
    bool requires_time = false;
    
    std::mt19937 generator;
    std::uniform_real_distribution<real> phase_distribution;
};
