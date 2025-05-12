#ifndef SYSTEMSOLVER_H
#define SYSTEMSOLVER_H

#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/Sparse>
#include <eigen3/Eigen/src/Core/Matrix.h>
#include <iostream>
#include <vector>
#include <umfpack.h>

using Eigen::MatrixXd;
using Eigen::VectorXd;

class SystemSolver
{
    public:
        //resolve sistema inicial
        SystemSolver(const std::vector<int>& basics_variables, const Eigen::SparseMatrix<double>& A);

        
        std::vector<std::pair<int, VectorXd>> E;
        Eigen::SparseMatrix<double> B0;

        //UMFPACK
        double *null;
        void *Symbolic, *Numeric;
        void refactor(const std::vector<int>& basics_variables, const Eigen::SparseMatrix<double>& A);
        VectorXd solve_initial(const VectorXd& RHS);
        VectorXd solve_price(const VectorXd& c_b);
        VectorXd solve_direction(const VectorXd& a);
        
};

#endif