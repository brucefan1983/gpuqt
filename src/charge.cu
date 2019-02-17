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


#include "charge.h"
#include <iostream>
#include <fstream>
#include <limits.h>
const double cutoff = 5.0;


void Charge::create_random_numbers
(std::mt19937& generator, int max_value, int total_number, int* random_numbers)
{
    int *permuted_numbers = new int[max_value];
    for(int i = 0; i < max_value; ++i)
    {
        permuted_numbers[i] = i;
    }
    std::uniform_int_distribution<int> rand_int(0, INT_MAX);
    for(int i = 0; i < max_value; ++i)
    {
        int j = rand_int(generator) % (max_value - i) + i;
        int temp = permuted_numbers[i];
        permuted_numbers[i] = permuted_numbers[j];
        permuted_numbers[j] = temp;
    }
    for (int i = 0; i < total_number; ++i)
    {
        random_numbers[i] = permuted_numbers[i];
    }
    delete[] permuted_numbers;
}


void Charge::find_potentials
(
    int number_of_atoms, double box_length[3], int pbc[3],
    std::vector<double>& x,  std::vector<double>& y,  std::vector<double>& z, 
    double* potential
)
{
    double cutoff_square = cutoff * cutoff * xi * xi;
    double xi_factor = -0.5 / (xi * xi);
    double box_length_half[3];
    for (int d = 0; d < 3; ++d) box_length_half[d] = box_length[d] * 0.5;
    for (int n = 0; n < number_of_atoms; ++n) potential[n] = 0.0;
    for (int i = 0; i < Ni; ++i)
    {
        int n1 = impurity_indices[i];
        double x1 = x[n1];
        double y1 = y[n1];
        double z1 = z[n1];
        for (int n2 = 0; n2 < number_of_atoms; ++n2)
        {
            double r12[3];
            r12[0] = x[n2] - x1;
            r12[1] = y[n2] - y1;
            r12[2] = z[n2] - z1;
            double d12_square = 0.0;
            for (int d = 0; d < 3; ++d)
            {
                r12[d] = fabs(r12[d]);
                if (pbc[d] == 1 && r12[d] > box_length_half[d])
                {
                    r12[d] = box_length[d] - r12[d];
                }
                d12_square += r12[d] * r12[d];
            }
            if (d12_square > cutoff_square) continue;
            potential[n2] += impurity_strength[i] * exp(d12_square * xi_factor);
        }
    }
}


void Charge::add_impurities
(
    std::mt19937& generator, int number_of_atoms, double box_length[3],
    int pbc[3],  std::vector<double>& x,  std::vector<double>& y, 
     std::vector<double>& z, double* potential
)
{
    impurity_indices = new int[Ni];
    impurity_strength = new double[Ni];
    create_random_numbers(generator, number_of_atoms, Ni, impurity_indices);
    double W2 = xi * 0.5;
    std::uniform_real_distribution<double> strength(-W2, W2);
    for (int i = 0; i < Ni; ++i)
    {
        impurity_strength[i] = strength(generator);
    }
    find_potentials(number_of_atoms, box_length, pbc, x, y, z, potential);
    delete[] impurity_indices;
    delete[] impurity_strength;
}


