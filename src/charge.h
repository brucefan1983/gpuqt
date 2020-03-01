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


class Charge
{
public:
    void add_impurities
    (
        std::mt19937&, int, double*, int*,  std::vector<double>&,
        std::vector<double>&,  std::vector<double>&, double*
    );
    bool has = false;
    int Ni;       // number of impurities
    double W;     // impurity strength
    double xi;    // impurity range
private:
    int Nx, Ny, Nz, Nxyz; // number of cells
    double rc;            // cutoff distance for impurity potential
    double rc2;           // cutoff square
    std::vector<int> cell_count;
    std::vector<int> cell_count_sum;
    std::vector<int> cell_contents;
    std::vector<int> impurity_indices;
    std::vector<double> impurity_strength;
    void find_impurity_indices(std::mt19937&, int);
    void find_impurity_strength(std::mt19937&);
    void find_potentials
    (
        int, double*, int*, std::vector<double>&, std::vector<double>&,
        std::vector<double>&, double*
    );
    int find_cell_id(double, double, double, double);
    void find_cell_id(double, double, double, double, int&, int&, int&, int&);
    void find_cell_numbers(int*, double*);
    void find_cell_contents
    (
        int, int*, double*, std::vector<double>&, 
        std::vector<double>&, std::vector<double>&
    );
    int find_neighbor_cell(int, int, int, int, int, int, int);
};


