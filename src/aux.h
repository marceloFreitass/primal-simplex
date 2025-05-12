#include <vector>
#include <iostream>
#include "Simplex.h"
#include "SystemSolver.h"
#include "mpsReader.h"

#define EPSILON_2 0.00000001
#define EPSILON_3 0.000001

void print_vector_int(const std::vector<int>& v);
void solve(mpsReader& Data, Simplex& simplex, SystemSolver& system_solver);
void solve_first_phase(mpsReader& Data, Simplex& simplex, SystemSolver& system_solver);
void update_basis_info(Simplex& simplex, VectorXd& d, int entering_variable, 
    int leaving_variable);
void update_first_phase(Simplex& simplex, VectorXd& d, int entering_variable, 
    int leaving_variable);
