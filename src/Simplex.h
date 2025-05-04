#ifndef SIMPLEX_H
#define SIMPLEX_H

#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/Sparse>
#include <eigen3/Eigen/src/Core/Matrix.h>
#include <iostream>
#include <vector>
#include <numeric>
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
        VectorXd c_b;
        VectorXd x;
        
        std::vector<int> basics;
        //mapeia os indices das variaveis basicas, se for uma variavel nao basica, = -1
        std::vector<int> basics_idx;

        std::vector<int> non_basics;
        //mapeia os indices das variaveis nao_basica, se for uma variavel basica, = -1
        std::vector<int> non_basics_idx;
        int n;
        int m;

        std::vector<int> violate_bound; //se v[i] = -1, x_i viola LB, se v[i] = 1 viola UB, 0 viavel

        double check_initial_infeasible(const VectorXd& init_x_b, const VectorXd& init_x_n);
        //variable and reduced cost
        std::pair<int,double> entering_variable(const VectorXd& y);
        //variavel, bool = 1 sse  a variavel que sai estava na base (a variavel nao basica so troca de bound) e valor de t;
        std::tuple<int, bool, double> leaving_variable(const VectorXd& d, const std::pair<int,double>& entering_variable);
        inline void update_sol(const VectorXd& d, const double t, const int non_basic_variable)
        {
            for(int i = 0; i < m; i++)
            {
                x(basics[i]) += t * d(i);
            }
            x(non_basic_variable) -= t;
        }
};

#endif