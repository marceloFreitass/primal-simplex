#include "aux.h"

void print_vector_int(const std::vector<int>& v)
{
    std::cout << "VEC: ";
    for(size_t i = 0; i < v.size(); i++)
    {
        std::cout << v[i] << " ";
    }
    std::cout << std::endl;
}

void solve(const mpsReader& Data)
{
    Simplex solver(Data.A, Data.lb, Data.ub, Data.c);
    SystemSolver SS(solver.basics, solver.A);
    VectorXd y;
    VectorXd d;
    double obj = 0;

    while(true)
    {
        y = SS.solve_price(solver.c_b);
        std::pair<int, double> enter_var = solver.entering_variable(y);
        if(enter_var.first == solver.n) // nenhum custo reduzido melhora
            break;
        
        d = SS.solve_direction(solver.A.col(enter_var.first));
            
        getchar();
    }
}