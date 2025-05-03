#include "Simplex.h"

Simplex::Simplex(const MatrixXd& A, const VectorXd& lb, const VectorXd& ub, const VectorXd& c)
{
    std::cout << "Matriz original: \n";
    std::cout << A << std::endl;
    this->A = A.sparseView();
    std::cout << "Esparsa: \n";
    std::cout << this->A << std::endl;
    this->lb = lb;
    this->ub = ub;
    this->c = c;

    n = A.cols(); //ja ta com de folga, numero real de variaives = n - m
    m = A.rows();
    basics.reserve(m);
    non_basics.reserve(n - m);
    c_b = VectorXd::Zero(m);
    std::cout << c_b << std::endl;
    std::cout << "N: " << n << std::endl;
    std::cout << "M: " << m << std::endl;
    x = VectorXd::Zero(n); // TODO: CONSEGUIR SOLUCAO INICIAL DE FORMA ALGORITMICA
    for(int i = n - m; i < n; i++)
    {
        std::cout << "ind: " << i - (n - m) << "\n";
        basics.push_back(i);
        c_b(i - (n - m)) = c(basics[i - (n - m)]);
        // aux++;
        // std::cout << c_b << std::endl;
    }
    for(int i = 0; i < n - m; i++)
    {
        non_basics.push_back(i);
    }

    print_vector_int(basics);
    print_vector_int(non_basics);

    std::cout << c_b << std::endl;


    

}

std::pair<int,double> Simplex::entering_variable(const VectorXd& y)
{
    int entering_variable = n;
    double var_reduced_cost = 0;
    for(size_t i = 0; i < non_basics.size(); i++)
    {   
        int var = non_basics[i];
        double reduced_cost = c(var) - y.transpose()  * A.col(var);
        
        // std::cout << "I: " << var << std::endl;
        // std::cout << "C: " << c(var) << std::endl;
        // std::cout << "RC: " << reduced_cost << std::endl;
        if(reduced_cost > EPSILON_1 && ub(var) - x(var) > EPSILON_1)  //custo reduzido positivo e pode aumentar
        {
            if(var < entering_variable)
            {
                entering_variable = var;
                var_reduced_cost = reduced_cost;
            }
        }
        else if(reduced_cost < -EPSILON_1 && lb(var) - x(var) < -EPSILON_1)  //custo reduzido negativo e pode diminuir
        {
            if(var < entering_variable)
            {
                entering_variable = var;
                var_reduced_cost = reduced_cost;
            }
        }
    }
    std::cout << "VAI ENTRAR: " << entering_variable << std::endl;
    std::cout << "CUSTO REDUZIDO: " << var_reduced_cost << std::endl;
    return {entering_variable, var_reduced_cost};
    
}

std::pair<int, bool> Simplex::leaving_variable(const VectorXd& d, const std::pair<int,double>& entering_variable)
{
    bool basic = false;
    int leaving_variable = entering_variable.first;


    
}
