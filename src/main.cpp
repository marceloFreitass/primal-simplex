#include <iostream>

#include "mpsReader.h"
#include "aux.h"

int main(int argc, char** argv)
{
    if(argc != 2)
    {
        std::cout << "Usage: ./exe file.mps";
    }

    mpsReader leitor(argv[1]);

    Simplex solver(leitor.A, leitor.lb, leitor.ub, leitor.c);
    SystemSolver SS(solver.basics, solver.A);
    std::cout << leitor.lb.transpose() << std::endl;
    std::cout << leitor.ub.transpose() << std::endl;
    //os problemas sao de minimizacao
    leitor.c = -leitor.c;

    solve_first_phase(leitor, solver, SS);
    solve(leitor, solver, SS);

    
    return 0;

}