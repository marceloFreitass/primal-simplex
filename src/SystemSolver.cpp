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
    // std::cout << "B: \n";

    // std::cout << B0 << std::endl;
    std::cout << "Tamanho de E: " << E.size() << std::endl;
    //DECOMPOSICAO LU DE B0
    null = (double *)NULL;

    (void)umfpack_di_symbolic(m, m, B0.outerIndexPtr(), B0.innerIndexPtr(), B0.valuePtr(), &Symbolic, null, null);
    (void)umfpack_di_numeric(B0.outerIndexPtr(), B0.innerIndexPtr(), B0.valuePtr(), Symbolic, &Numeric, null, null);
}

VectorXd SystemSolver::solve_initial(const VectorXd& RHS)
{
    int m = RHS.size();
    VectorXd x_b = VectorXd::Zero(m);
    
    (void)umfpack_di_solve(UMFPACK_A, B0.outerIndexPtr(), B0.innerIndexPtr(), B0.valuePtr(), x_b.data(), RHS.data(), Numeric, null, null);
    
    return x_b;
}

VectorXd SystemSolver::solve_price(const VectorXd& c_b)
{
    VectorXd y = c_b;
    std::cout << "C_b: " << c_b.transpose() << std::endl;
    int k = E.size();
    int m = c_b.size();
    //TODO : CASO EM QUE É 0
    for(int i = k - 1; i >= 0; i--)
    {
        int p = E[i].first;
        for(int j = 0; j < p; j++)
        {
            y(p) -= (y(j) * E[i].second(j));
        }
        for(int j = p + 1; j < m; j++)
        {
            y(p) -= (y(j) * E[i].second(j));

        }
        y(p) /= E[i].second(p);
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

    int k = E.size();
    for(int i = 0; i < k; i++)
    {
        int p = E[i].first;
        d(p) /= E[i].second(p);
        // std::cout << d(p) << std::endl;
        for(int j = 0; j < p; j++)
        {
            d(j) -= d(p) * E[i].second(j); 
        }
        for(int j = p + 1; j < m; j++)
        {
            d(j) -= d(p) * E[i].second(j); 
        }
    }
    return d;
}