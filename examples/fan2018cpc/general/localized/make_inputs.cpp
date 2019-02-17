#include <stdlib.h>
#include <vector>
#include <stdio.h>
#include <time.h>
#include <random>
#include <chrono>


// some parameters
#define Nx 4000000               
#define Ny 1              
#define N  (Nx * Ny) // number of sites
#define MN 2         // square lattice
#define pbc_x 1      // periodic
#define pbc_y 0      // periodic
#define W 4.0        // Anderson disorder strength


int find_index(int nx, int ny)
{
    if (nx < 0) nx += Nx;
    if (nx >= Nx) nx -= Nx;
    if (ny < 0) ny += Ny;
    if (ny >= Ny) ny -= Ny;		
    int index = nx * Ny + ny; 
    return index;
}


void generate_square_lattice(std::vector<int>& NN, std::vector<int>& NL, std::vector<double>& x)
{
    for (int i = 0; i < Nx; ++i)
    {
        for (int j = 0; j < Ny; ++j)
        {
            int index = find_index(i, j);	
            x[index] = i;
            NN[index] = 0;
            if (pbc_x || i > 0) // left atom
            {
                NL[MN * index + NN[index]] = find_index(i - 1, j);
                NN[index]++;
            }
            if (pbc_x || i < Nx - 1) // right atom
            {
                NL[MN * index + NN[index]] = find_index(i + 1, j);
                NN[index]++;
            }
            if (pbc_y || j < Ny - 1) // atom above
            {
                NL[MN * index + NN[index]] = find_index(i, j + 1);
                NN[index]++;
            }
            if (pbc_y || j > 0) // atom below
            {
                NL[MN * index + NN[index]] = find_index(i, j - 1);
                NN[index]++;
            }			
        }
    }
}


void generate_anderson_disorder(std::vector<double>& potential)
{
    std::mt19937 generator(std::chrono::system_clock::now().time_since_epoch().count());
    std::uniform_real_distribution<double> real_distribution(-W * 0.5, W * 0.5);
    for (int i = 0; i < N; ++i)
    {
        potential[i] = real_distribution(generator);
    }
}


int main(int argc, char* argv[])
{
    srand(time(0));
	
    std::vector<int> NN, NL;	
    std::vector<double> x, potential; 
    NN.resize(N);
    NL.resize(N * MN);
    x.resize(N);
    potential.resize(N);
	
    generate_square_lattice(NN, NL, x);
    generate_anderson_disorder(potential);
	
    // open files for output
    FILE *fid_energy    = fopen("energy.in", "w");
    FILE *fid_time_step = fopen("time_step.in", "w");
    FILE *fid_para      = fopen("para.in", "w");
    FILE *fid_neighbor  = fopen("neighbor.in",  "w");
    FILE *fid_position  = fopen("position.in",  "w");
    FILE *fid_potential = fopen("potential.in", "w");

    // energy.in
    fprintf(fid_energy, "%d\n", 401);
    for (int n = -200; n <= 200; n++)
    {
        fprintf(fid_energy, "%f\n", n * 0.02);
    }

    // time_step.in (increasingly large time steps)
    double base = 0.1;
    fprintf(fid_time_step, "%d\n", 30);
    for (int n1 = 0; n1 < 3; n1++)
    {
        for (int n2 = 1; n2 <= 10; n2++)
        {
            fprintf(fid_time_step, "%f\n", n2 * base);
        }
        base *= 10;
    }
    
    // para.in
    fprintf(fid_para, "number_of_random_vectors 2\n");
    fprintf(fid_para, "number_of_moments        500\n");
    fprintf(fid_para, "energy_max               4.1\n");
    fprintf(fid_para, "calculate_msd\n");

    // neighbor.in
    fprintf(fid_neighbor, "%10d%10d\n", N, MN);
    for (int n = 0; n < N; ++n)
    {
        fprintf(fid_neighbor, "%10d", NN[n]); 
        for (int m = 0; m < NN[n]; ++m)
        {
            fprintf(fid_neighbor, "%10d", NL[n * MN + m]);
        }
        fprintf(fid_neighbor, "\n");
    }

    // position.in
    fprintf(fid_position, "%15.6e%15.6e\n", double(Nx), double(Nx) * double(Ny));
    for (int n = 0; n < N; ++n)
    {
        fprintf(fid_position, "%15.6e\n", x[n]);
    }

    // potential.in
    for (int n = 0; n < N; ++n)
    {
        fprintf(fid_potential, "%15.6e\n", potential[n]);
    }
    
    // close the files
    fclose(fid_energy);
    fclose(fid_time_step);
    fclose(fid_para);	
    fclose(fid_neighbor);
    fclose(fid_position);
    fclose(fid_potential);
    	
    //system("PAUSE"); for DEV C++
    return 0;
}


