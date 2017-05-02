#pragma once

#include "common.h"
#include "vector.h"
#include "hamiltonian.h"

class Model;

void find_dos(Model &para, Hamiltonian& H, Vector& random_state);
void find_vac(Model &para, Hamiltonian& H, Vector& random_state);
void find_msd(Model &para, Hamiltonian& H, Vector& random_state);


