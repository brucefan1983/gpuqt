#pragma once
#include "common.h"

class Model;


class Vector
{
 
public:
    Vector(int n, Model& parameters);
    Vector(Model& parameters);
    Vector(Vector& original);
    Vector(real* original_real, real* original_imag, Model& parameters, bool device=true);
	
    ~Vector();
    void add(Vector& other, real coeff=1.0);
    void copy(Vector& other);
    void copy_from_host(real* other_real, real* other_imag);
    void copy_to_host(real* target_real, real* target_imag);	
    void swap(Vector& other);
    void inner_product_1(Vector& other, Vector& target, int offset);
    void inner_product_2(Vector& target);	
    
    real* real_part;
    real* imag_part;
    
private:
    void initialize_parameters();
    int n;
    size_t array_size;
    size_t grid_size;
    Model& model;
    
};

