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
#include <random>
#include <iostream>
#include <fstream>




void Model::add_anderson_disorder()
{
    potential = new real[number_of_atoms];
    real W2 = anderson_disorder_strength * 0.5;
    std::uniform_real_distribution<real> on_site_potential(-W2, W2);
    for (int n = 0; n < number_of_atoms; ++n)
    {
        if (has_anderson_disorder)
            potential[n] = on_site_potential(generator);
        else
            potential[n] = 0.0;
    }
}




void Model::create_random_numbers
(int max_value, int total_number, int* random_numbers)
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




void Model::specify_vacancies
(int *is_vacancy, int number_of_atoms_pristine)
{
    int *vacancy_indices = new int[number_of_vacancies];
    create_random_numbers
    (number_of_atoms_pristine, number_of_vacancies, vacancy_indices);

    for (int n = 0; n < number_of_atoms_pristine; ++n)
    {
        is_vacancy[n] = 0;
    }
    for (int n = 0; n < number_of_vacancies; ++n)
    {
        is_vacancy[vacancy_indices[n]] = 1;
    }
    delete[] vacancy_indices;
}




void Model::find_new_atom_index
(int *is_vacancy, int *new_atom_index, int number_of_atoms_pristine)
{
    int count = 0;
    for (int n = 0; n < number_of_atoms_pristine; ++n)
    {
        if (is_vacancy[n] == 0)
        {
            new_atom_index[n] = count;
            ++count;
        }
    }
}




void Model::add_vacancies()
{
    // copy some data
    int *neighbor_number_pristine = new int[number_of_atoms];
    int *neighbor_list_pristine = new int[number_of_pairs];
    real *hopping_real_pristine = new real[number_of_pairs];
    real *hopping_imag_pristine = new real[number_of_pairs];
    real *xx_pristine = new real[number_of_pairs];

    for (int n = 0; n < number_of_atoms; ++n)
    {
        neighbor_number_pristine[n] = neighbor_number[n];
    }   
    for (int m = 0; m < number_of_pairs; ++m)
    {
        neighbor_list_pristine[m] = neighbor_list[m];
        hopping_real_pristine[m] = hopping_real[m];
        hopping_imag_pristine[m] = hopping_imag[m];
        xx_pristine[m] = xx[m];
    }

    // change parameters
    int number_of_atoms_pristine = number_of_atoms;
    number_of_atoms = number_of_atoms_pristine - number_of_vacancies;
    number_of_pairs = number_of_atoms * max_neighbor;

    // delete old memory
    delete[] neighbor_number;
    delete[] neighbor_list;
    delete[] hopping_real;
    delete[] hopping_imag;
    delete[] xx;

    // allocate new memory
    neighbor_number = new int[number_of_atoms];
    neighbor_list = new int[number_of_pairs];
    hopping_real = new real[number_of_pairs];
    hopping_imag = new real[number_of_pairs];
    xx = new real[number_of_pairs];

    // specify the distribution of the vacancies
    int *is_vacancy = new int[number_of_atoms_pristine];
    specify_vacancies(is_vacancy, number_of_atoms_pristine);

    // find the new indices of the atoms
    int *new_atom_index = new int[number_of_atoms_pristine];
    find_new_atom_index(is_vacancy, new_atom_index, number_of_atoms_pristine);

    // get the new neighbor structure and related data
    int count_atom = 0;
    for (int n = 0; n < number_of_atoms_pristine; ++n)
    {
        if (is_vacancy[n] == 0)
        {
            int count_neighbor = 0;
            for (int m = 0; m < neighbor_number_pristine[n]; ++m)
            {
                int index_old = n * max_neighbor + m;
                int k = neighbor_list_pristine[index_old];
                if (is_vacancy[k] == 0)
                {
                    int index_new = count_atom * max_neighbor + count_neighbor;
                    neighbor_list[index_new] = new_atom_index[k];
                    hopping_real[index_new] = hopping_real_pristine[index_old];
                    hopping_imag[index_new] = hopping_imag_pristine[index_old];
                    xx[index_new] = xx_pristine[index_old];
                    ++count_neighbor;
                }
            }
            neighbor_number[count_atom] = count_neighbor;
            ++count_atom;
        }
    }

    // free memory
    delete[] neighbor_number_pristine;
    delete[] neighbor_list_pristine;
    delete[] hopping_real_pristine;
    delete[] hopping_imag_pristine;
    delete[] xx_pristine;
    delete[] is_vacancy;  
    delete[] new_atom_index;
}




static int find_index
(int nx, int ny, int nz, int Nx, int Ny, int Nz, int m, int N_orbital)
{
    if (nx < 0) nx += Nx;
    if (nx >= Nx) nx -= Nx;
    if (ny < 0) ny += Ny;
    if (ny >= Ny) ny -= Ny;
    if (nz < 0) nz += Nz;
    if (nz >= Nz) nz -= Nz;
    return ((nx * Ny + ny) * Nz + nz) * N_orbital + m;
}




void Model::initialize_lattice_model()
{
    std::string filename = input_dir + "/lattice.in";
    print_started_reading(filename);
    std::ifstream input(filename);

    if (!input.is_open())
    {
        std::cout <<"Could not open " + filename << std::endl;
        exit(1);
    }
 
    int N_orbital;
    int transport_direction;
    int N_cell[3];
    int pbc[3];
    real box[3];
    real lattice_constant[3];

    input >> N_cell[0] >> N_cell[1] >> N_cell[2];
    std::cout << "number of cells  = " 
         << N_cell[0] << " " << N_cell[1] << " " << N_cell[2] << std::endl;

    input >> pbc[0] >> pbc[1] >> pbc[2] >> transport_direction;
    std::cout << "pbc = " << pbc[0] << " " << pbc[1] << " " << pbc[2] 
              << std::endl;
    std::cout << "transport direction = " << transport_direction << std::endl;

    if (pbc[transport_direction] != 1)
    {
        std::cout << "Error: transport direction must be periodic" << std::endl;
        exit(1);
    }

    input >> lattice_constant[0] >> lattice_constant[1] >> lattice_constant[2];
    std::cout << "lattice constant = " 
         << lattice_constant[0] << " "
         << lattice_constant[1] << " "
         << lattice_constant[2] << " "
         << std::endl;
    for (int d = 0; d < 3; ++d)
        box[d] = lattice_constant[d] * N_cell[d];
    volume = box[0] * box[1] * box[2];
    std::cout << "box = " << box[0] << " " << box[1] << " " << box[2] << " "
              << std::endl;

    input >> N_orbital >> max_neighbor;
    std::cout << "nnumber of orbitals per cell = " << N_orbital << std::endl;
    std::cout << "maximum number of hoppings per orbital = " << max_neighbor
              << std::endl;
    number_of_atoms = N_orbital * N_cell[0] * N_cell[1] * N_cell[2];
    std::cout << "number_of_atoms = " << number_of_atoms << std::endl;

    number_of_pairs = number_of_atoms * max_neighbor;
    neighbor_number = new int[number_of_atoms];
    neighbor_list = new int [number_of_pairs];
    hopping_real = new real[number_of_pairs];
    hopping_imag = new real[number_of_pairs];
    xx = new real[number_of_pairs];

    std::vector<real> x_cell;
    x_cell.resize(N_orbital);
    int number_of_hoppings_per_cell = N_orbital * max_neighbor;
    std::vector<std::vector<int>> hopping_data;
    hopping_data.assign(6, std::vector<int>(number_of_hoppings_per_cell, 0));

    std::cout << std::endl << "orbital\tx" << std::endl;
    for (int n = 0; n < N_orbital; ++n)
    {
        input >> x_cell[n];
        std::cout << n << "\t" << x_cell[n] << std::endl;
    }

    std::vector<int> number_of_hoppings;
    number_of_hoppings.resize(N_orbital);
    for (int m = 0; m < N_orbital; m++)
    {
        input >> number_of_hoppings[m];
        std::cout << std::endl << "number_of_hoppings for orbital " << m
                  << " = " << number_of_hoppings[m] << std::endl;

        for (int n = 0; n < number_of_hoppings[m]; ++n)
        {
            int nx, ny, nz, m_neighbor;
            real hopping_real, hopping_imag;
            input >> nx >> ny >> nz >> m_neighbor >> hopping_real
                  >> hopping_imag;

            hopping_data[0][m*max_neighbor+n] = nx;
            hopping_data[1][m*max_neighbor+n] = ny;
            hopping_data[2][m*max_neighbor+n] = nz;
            hopping_data[3][m*max_neighbor+n] = m_neighbor;
            hopping_data[4][m*max_neighbor+n] = hopping_real;
            hopping_data[5][m*max_neighbor+n] = hopping_imag;

            std::cout << "H(0,0,0," << m << "; " 
                 << nx << "," << ny << "," << nz << "," << m_neighbor << ") = "
                 << hopping_real << " + i " << hopping_imag << std::endl;
        }
    }

    for (int nx1 = 0; nx1 < N_cell[0]; ++nx1)
    {
        for (int ny1 = 0; ny1 < N_cell[1]; ++ny1)
        {  
            for (int nz1 = 0; nz1 < N_cell[2]; ++nz1)
            {
                for (int m = 0; m < N_orbital; ++m)
                {
                    int n1 = find_index
                    (
                        nx1, ny1, nz1, N_cell[0], N_cell[1], N_cell[2], 
                        m, N_orbital
                    );

                    int count = 0;
                    for (int i = 0; i < number_of_hoppings[m]; ++i)
                    {
                        int neighbor_index = n1 * max_neighbor + count;
                        int k = m*max_neighbor+i;

                        int nx2 = hopping_data[0][k] + nx1;
                        int ny2 = hopping_data[1][k] + ny1;
                        int nz2 = hopping_data[2][k] + nz1;
                        bool skip_x = !pbc[0] && (nx2 < 0 || nx2 >= N_cell[0]);
                        bool skip_y = !pbc[1] && (ny2 < 0 || ny2 >= N_cell[1]);
                        bool skip_z = !pbc[2] && (nz2 < 0 || nz2 >= N_cell[2]);
                        if (skip_x || skip_y || skip_z) continue;

                        neighbor_list[neighbor_index] = find_index
                        (
                            nx2, ny2, nz2, N_cell[0], N_cell[1], N_cell[2], 
                            hopping_data[3][k], N_orbital
                        );

                        real x12 = lattice_constant[transport_direction]
                                 * hopping_data[transport_direction][k];
                        x12 += x_cell[hopping_data[3][k]] - x_cell[m];
                        xx[neighbor_index] = x12;

                        hopping_real[neighbor_index] = hopping_data[4][k];
                        hopping_imag[neighbor_index] = hopping_data[5][k];

                        ++count;
                    } 
                    neighbor_number[n1] = count;
                }
            }
        }
    }

    if (has_vacancy_disorder && number_of_vacancies > 0)
    {
        add_vacancies();
    }
    add_anderson_disorder();
    print_finished_reading(filename);
}



