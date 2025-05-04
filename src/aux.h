#include <vector>
#include <iostream>
#include "Simplex.h"
#include "SystemSolver.h"
#include "mpsReader.h"

#define EPSILON_2 0.00000001
#define EPSILON_3 0.000001

void print_vector_int(const std::vector<int>& v);
void solve(const mpsReader& Data);
void update_basis_info(Simplex& simplex, const VectorXd& d, const int entering_variable, 
    const int leaving_variable);
void update_first_phase(Simplex& simplex, const VectorXd& d, const int entering_variable, 
    const int leaving_variable);
void solve_first_phase(const mpsReader& Data);