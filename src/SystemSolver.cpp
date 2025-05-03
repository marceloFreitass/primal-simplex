#include "SystemSolver.h"

//checar se a matriz vira um atributo
SystemSolver::SystemSolver(const std::vector<int>& basics_variables, const Eigen::SparseMatrix<double>& A)
{

    int m = basics_variables.size();
    B0 = Eigen::SparseMatrix<double> (m, m);
    for(int i = 0; i < m; i++)
    {
        B0.col(i) = A.col(basics_variables[i]);
    }
    std::cout << "B: \n";

    std::cout << B0 << std::endl;
    std::cout << "Tamanho de E: " << E.size() << std::endl;
    //DECOMPOSICAO LU DE B0
    null = (double *)NULL;

    (void)umfpack_di_symbolic(m, m, B0.outerIndexPtr(), B0.innerIndexPtr(), B0.valuePtr(), &Symbolic, null, null);
    (void)umfpack_di_numeric(B0.outerIndexPtr(), B0.innerIndexPtr(), B0.valuePtr(), Symbolic, &Numeric, null, null);
}

VectorXd SystemSolver::solve_price(const VectorXd& c_b)
{
    VectorXd y = c_b;
    int k = E.size();
    int m = c_b.size();
    //TODO: if se 0
    //OLD: errado (transposto é diferente)
    // for(int i = k - 1; i >= 0; i--)
    // {

    //     int p = E[i].first;
    //     y(p) /= E[i].second(p);
    //     for(int j = 0; i < m; i++)
    //     {
    //         if(j == p)
    //         {
    //             continue;
    //         }
    //         y(i) -= y(p) * E[i].second(j);
    //     }
    // }

    for(int i = k - 1; i >= 0; i--)
    {
        int p = E[i].first;
        y(p) /= E[i].second(p);
        for(int j = 0; j < p; j++)
        {
            y(p) -= y(j)/E[i].second(j);
        }
        for(int j = p; p < m; j++)
        {
            y(p) -= y(j)/E[i].second(j);
        }
    }
    VectorXd result = VectorXd::Zero(m);
    (void)umfpack_di_solve(UMFPACK_At, B0.outerIndexPtr(), B0.innerIndexPtr(), B0.valuePtr(), result.data(), y.data(), Numeric, null, null);
    return result;
}

VectorXd SystemSolver::solve_direction(const VectorXd& a)
{
    int m = a.size();
    VectorXd d = VectorXd::Zero(m);
    (void)umfpack_di_solve(UMFPACK_A, B0.outerIndexPtr(), B0.innerIndexPtr(), B0.valuePtr(), d.data(), a.data(), Numeric, null, null);
    std::cout << d << std::endl;


    int k = E.size();
    // int m = a.size().
    for(int i = 0; i < k; i++)
    {
        int p = E[i].first;
        d(p) /= E[i].second(p);
        for(int j = 0; j < p; j++)
        {
            d(j) -= d(p) * E[i].second(j); 
        }
        for(int j = p; j < m; j++)
        {
            d(j) -= d(p) * E[i].second(j); 
        }
    }
    return d;
}