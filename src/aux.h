#include <vector>
#include <iostream>
#include "Simplex.h"
#include "SystemSolver.h"
#include "mpsReader.h"

#define EPSILON_2 0.00000001
#define EPSILON_3 0.000001

void print_vector_int(const std::vector<int>& v);
void solve(const mpsReader& Data);