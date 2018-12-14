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

 
#include "vector.h"
#include "model.h"


/*
	Gets number of elements from model and sets array and grid sizes accordingly.
	Allocates memory on the device
*/
void Vector::initialize_parameters()
{
    n = model.number_of_atoms;
    array_size = n * sizeof(real);
    grid_size = (n-1) / BLOCK_SIZE + 1;
    cudaMalloc((void**)&real_part, array_size);
    cudaMalloc((void**)&imag_part, array_size);	
}


/*
	Kernel for setting all elements of a state to zero (both real and imaginary parts)
*/
__global__ void gpu_set_zero(int number_of_elements, real *g_state_real, real *g_state_imag)
{
    int n = blockIdx.x * blockDim.x + threadIdx.x;
    if (n < number_of_elements)
    {
        g_state_real[n] = 0; 
        g_state_imag[n] = 0;  
    }
}


/*
	Constructor for an empty vector. 
	Takes length from model and sets all elements to zero
*/
Vector::Vector(Model& model) : model(model)
{
    initialize_parameters();	
    gpu_set_zero<<<grid_size, BLOCK_SIZE>>>(n, real_part, imag_part);
}



/*
	Constructor for a vector of arbitrary length. 
	Does not initialize data. 
*/
Vector::Vector(int n, Model& model) : model(model)
{
    this->n = n;
    array_size = n * sizeof(real);
    grid_size = (model.number_of_atoms-1) / BLOCK_SIZE + 1;
    cudaMalloc((void**)&real_part, array_size);
    cudaMalloc((void**)&imag_part, array_size);
}



/*
	Kernel for copying states on the gpu
*/
__global__ void gpu_copy_state
(int number_of_atoms, real *g_state_in_real, real *g_state_in_imag, real *g_state_out_real, real *g_state_out_imag)
{
    int n = blockIdx.x * blockDim.x + threadIdx.x;
    if (n < number_of_atoms)
    {
        g_state_out_real[n] = g_state_in_real[n]; 
        g_state_out_imag[n] = g_state_in_imag[n];  
    }
}


/*
	Constructor which creates a copy of *original*
*/
Vector::Vector(Vector& original) : model(original.model)
{
    initialize_parameters();
    gpu_copy_state<<<grid_size, BLOCK_SIZE>>>(n, original.real_part, original.imag_part, real_part, imag_part);
}




// Destructor
Vector::~Vector()
{
    cudaFree(real_part);
    cudaFree(imag_part);
}


// Add the "other" vector to the current vector. This is the kernel
__global__ void gpu_add_state
(int n, real *g_state_in_real, real *g_state_in_imag, real *g_state_out_real, real *g_state_out_imag)
{
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i < n)
    {
        g_state_out_real[i] += g_state_in_real[i]; 
        g_state_out_imag[i] += g_state_in_imag[i];  
    }
}


// Add the "other" vector to the current vector. This is a wrapper function
void Vector::add(Vector& other, real coeff)
{
    gpu_add_state<<<grid_size, BLOCK_SIZE>>>(n, other.real_part, other.imag_part, real_part, imag_part);
}


// Sets this vector to the same state as "other" vector
void Vector::copy(Vector& other)
{
    if (other.n == n)
    {
        gpu_copy_state<<<grid_size, BLOCK_SIZE>>>(n, other.real_part, other.imag_part, real_part, imag_part);	
    }
    else
    {
        std::cout << "Array sizes do not match in copy." << std::endl;	
    }
}


// Copy from a host vector to the current vector
void Vector::copy_from_host(real* other_real, real* other_imag)
{
    cudaMemcpy(real_part, other_real, array_size, cudaMemcpyHostToDevice);
    cudaMemcpy(imag_part, other_imag, array_size, cudaMemcpyHostToDevice);		
}


// Copy the current vector to a host vector
void Vector::copy_to_host(real* target_real, real* target_imag)
{
    cudaMemcpy(target_real, real_part, array_size, cudaMemcpyDeviceToHost);
    cudaMemcpy(target_imag, imag_part, array_size, cudaMemcpyDeviceToHost);
}


// Exchange the pointers of the current vector and the "other" vector
void Vector::swap(Vector& other)
{
    real* tmp_real = real_part;
    real* tmp_imag = imag_part;
    real_part = other.real_part,
    imag_part = other.imag_part;
    other.real_part = tmp_real;
    other.imag_part = tmp_imag;
}


/*
	Device function which performs sum reduction over warp
*/
__device__ void warp_reduce(volatile real *s, int t)
{
    s[t] += s[t + 32]; s[t] += s[t + 16]; s[t] += s[t + 8];
    s[t] += s[t + 4];  s[t] += s[t + 2];  s[t] += s[t + 1];
}


// The first step of calculating the inner products. This is the kernel
__global__ void gpu_find_inner_product_1
(
    int number_of_atoms,
    real *g_final_state_real, 
    real *g_final_state_imag, 
    real *g_random_state_real,
    real *g_random_state_imag,
    real *g_inner_product_real, 
    real *g_inner_product_imag,
    int   g_offset
)
{
    int tid = threadIdx.x;
    int n = blockIdx.x * blockDim.x + tid;
    int m;
    real a, b, c, d;
    __shared__ real s_data_real[BLOCK_SIZE];
    __shared__ real s_data_imag[BLOCK_SIZE];
    s_data_real[tid] = 0.0;
    s_data_imag[tid] = 0.0;
    
    if (n < number_of_atoms)
    {
        a = g_final_state_real[n];
        b = g_final_state_imag[n];
        c = g_random_state_real[n];
        d = g_random_state_imag[n];
        s_data_real[tid] = (a * c + b * d); 
        s_data_imag[tid] = (b * c - a * d);
    }
    __syncthreads();

    if (tid < 256) {m = tid + 256; s_data_real[tid] += s_data_real[m]; s_data_imag[tid] += s_data_imag[m];}
    __syncthreads();
    if (tid < 128) {m = tid + 128; s_data_real[tid] += s_data_real[m]; s_data_imag[tid] += s_data_imag[m];}
    __syncthreads();
    if (tid < 64)  {m = tid + 64;  s_data_real[tid] += s_data_real[m]; s_data_imag[tid] += s_data_imag[m];}
    __syncthreads();
    if (tid < 32)  {warp_reduce(s_data_real, tid); warp_reduce(s_data_imag, tid);}
    if (tid == 0) 
    {        
        g_inner_product_real[blockIdx.x + g_offset] = s_data_real[0];
        g_inner_product_imag[blockIdx.x + g_offset] = s_data_imag[0];
    }
}


// The first step of calculating the inner products. This is a wrapper function
void Vector::inner_product_1(Vector& other, Vector& target, int offset)
{
    gpu_find_inner_product_1<<<grid_size, 512>>>
    (
        model.number_of_atoms, real_part, imag_part, 
        other.real_part, other.imag_part, target.real_part, target.imag_part, 
        offset
    );
}


// The second step of calculating the inner products. This is the kernel
__global__ void gpu_find_inner_product_2
(
    int number_of_atoms,	
    real *g_inner_product_1_real, 
    real *g_inner_product_1_imag,
    real *g_inner_product_2_real, 
    real *g_inner_product_2_imag
)
{
    //<<<para.number_of_energy_points, BLOCK_SIZE)>>>
    int tid = threadIdx.x;
    int patch, n, m;

    __shared__ real s_data_real[BLOCK_SIZE];
    __shared__ real s_data_imag[BLOCK_SIZE];
    s_data_real[tid] = 0.0;
    s_data_imag[tid] = 0.0;
    int number_of_blocks  = (number_of_atoms - 1) / BLOCK_SIZE + 1;
    int number_of_patches = (number_of_blocks - 1) / BLOCK_SIZE + 1;

    for (patch = 0; patch < number_of_patches; ++patch)
    {
        n = tid + patch * BLOCK_SIZE;
        if (n < number_of_blocks)
        {
            m = blockIdx.x * number_of_blocks + n;
            s_data_real[tid] += g_inner_product_1_real[m]; 
            s_data_imag[tid] += g_inner_product_1_imag[m];
        }
    }
    __syncthreads();
  
    if (tid < 256) {m = tid + 256; s_data_real[tid] += s_data_real[m]; s_data_imag[tid] += s_data_imag[m];}
    __syncthreads();
    if (tid < 128) {m = tid + 128; s_data_real[tid] += s_data_real[m]; s_data_imag[tid] += s_data_imag[m];}
    __syncthreads();
    if (tid < 64) {m = tid + 64; s_data_real[tid] += s_data_real[m]; s_data_imag[tid] += s_data_imag[m];}
    __syncthreads();
    if (tid < 32) {warp_reduce(s_data_real, tid); warp_reduce(s_data_imag, tid);}
    if (tid == 0) 
    {        
        g_inner_product_2_real[blockIdx.x] = s_data_real[0];
        g_inner_product_2_imag[blockIdx.x] = s_data_imag[0];
    }    
}


// The second step of calculating the inner products. This is a wrapper function
void Vector::inner_product_2(Vector& target)
{
    gpu_find_inner_product_2<<<model.number_of_moments, 512>>>
    (
        model.number_of_atoms, real_part, imag_part, 
        target.real_part, target.imag_part
    );	
}



