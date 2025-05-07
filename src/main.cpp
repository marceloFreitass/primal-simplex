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
    // std::cout <<  "N ROWS EQ: " << leitor.n_rows_eq << std::endl;
    // std::cout << "N ROWS INEQ: " << leitor.n_rows_inq << std::endl;
    // std::cout << "N ROWS: " << leitor.n_rows << std::endl;
    // std::cout << "N VAR: " << leitor.n_cols << std::endl;
    // std::cout << leitor.b << std::endl;
    // std::cout << leitor.c << std::endl;
    std::cout << leitor.lb.transpose() << std::endl;
    std::cout << leitor.ub.transpose() << std::endl;
    //os problemas sao de minimizacao
    leitor.c = -leitor.c;

    solve_first_phase(leitor);

    return 0;
    std::cout << "WTF\n";

}