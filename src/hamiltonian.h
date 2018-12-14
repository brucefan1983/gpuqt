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
#include "common.h"
#include <chrono>




class Vector;
class Model;




class Hamiltonian
{
public:

    Hamiltonian(Model& para);
    ~Hamiltonian();

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
    //Model& model;
    int n;
    real energy_max;
};




