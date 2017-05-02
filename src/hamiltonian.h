#pragma once
#include "common.h"
#include <chrono>

class Vector;
class Model;

class Hamiltonian
{
 
public:
    Hamiltonian(Model& para);
    ~Hamiltonian();
    Hamiltonian
    (
        int* neighbor_number, int* neighbor_list, real* potential, 
        real* hopping_real, real* hopping_imag, Model& para, real* xx=NULL
    );
    void apply(Vector& input, Vector& output);
    void apply_commutator(Vector& input, Vector& output);
    void apply_current(Vector& input, Vector& output);
    void kernel_polynomial(Vector& state_0, Vector& state_1, Vector& state_2);
    void chebyshev_01(Vector& state_0, Vector& state_1, Vector& state, real bessel_0, real bessel_1, int direction);
    void chebyshev_2(Vector& state_0, Vector& state_1, Vector& state_2, Vector& state, real bessel_m, int label);
    void chebyshev_1x(Vector& input, Vector& output, real bessel_1);
    void chebyshev_2x
    (
        Vector& state_0, Vector& state_0x, Vector& state_1, Vector& state_1x, Vector& state_2, 
        Vector& state_2x, Vector& state, real bessel_m, int label
    );	

private:
    int* neighbor_number;
    int* neighbor_list;
    real* potential;
    real* hopping_real;
    real* hopping_imag;
    real* xx;
    size_t grid_size;
    Model& model;
    int n;
    real energy_max;
};




