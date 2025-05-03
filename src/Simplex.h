#ifndef Simplex_h
#define Simplex_h

#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/Sparse>
#include <eigen3/Eigen/src/Core/Matrix.h>
#include <iostream>
#include <vector>

#include "aux.h"
using Eigen::MatrixXd;
using Eigen::VectorXd;


#define EPSILON_1 0.00001 // 0

class Simplex
{
    public:
        Simplex(const MatrixXd& A, const VectorXd& lb, const VectorXd& ub, const VectorXd& c);


        Eigen::SparseMatrix<double> A;
        VectorXd lb;
        VectorXd ub;
        VectorXd c;
        //TODO: ficar tratando o c^TB (swapar junto com a base)
        VectorXd c_b;
        VectorXd x;
        std::vector<int> basics;
        std::vector<int> non_basics;
        int n;
        int m;

        //variable and reduced cost
        std::pair<int,double> entering_variable(const VectorXd& y);
        //variavel e bool = 1 sse  a variavel que sai estava na base (a variavel nao basica so troca de bound)
        std::pair<int, bool> leaving_variable(const VectorXd& d, const std::pair<int,double>& entering_variable);
        
};

#endif