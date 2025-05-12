#include "Simplex.h"

Simplex::Simplex(const MatrixXd& A, const VectorXd& lb, const VectorXd& ub, const VectorXd& c)
{

    this->A = A.sparseView();
    this->lb = lb;
    this->ub = ub;
    n = A.cols(); //ja ta com de folga, numero real de variaives = n - m
    m = A.rows();
    this->c = VectorXd::Zero(n);
    basics.reserve(m);
    non_basics.reserve(n - m);
    basics_idx = std::vector<int>(n, -1);
    non_basics_idx = std::vector<int>(n, -1);
    c_b = VectorXd::Zero(m);
    x = VectorXd::Zero(n);
    //forma B = I;
    for(int i = n - m; i < n; i++)
    {
        basics.push_back(i);
        c_b(i - (n - m)) = c(basics[i - (n - m)]);
        basics_idx[i] = (i - (n - m));
    }
    for(int i = 0; i < n - m; i++)
    {
        non_basics.push_back(i);
        non_basics_idx[i] = i;
    }
}

int Simplex::check_infeasible(const VectorXd& lb_, const VectorXd& ub_)
{
    int infeasible_variables = 0;
    double inf = std::numeric_limits<double>::infinity();
    c_b = VectorXd::Zero(m);
    c = VectorXd::Zero(n);
    ub = ub_;
    lb = lb_;
    for(int i = 0; i < m; i++)
    {
        if(lb_(basics[i]) - x(basics[i]) > epsilon)
        {
            c_b(i) = 1;
            infeasible_variables += 1;
            lb(basics[i]) = -inf;
            ub(basics[i]) = lb_(basics[i]);
        }
        else if(x(basics[i]) - ub_(basics[i]) > epsilon)
        {
            c_b(i) = -1;
            infeasible_variables += 1;
            ub(basics[i]) = inf;
            lb(basics[i]) = ub_(basics[i]);
        }
    }
    
    return infeasible_variables;
}

std::pair<int,double> Simplex::entering_variable(const VectorXd& y)
{
    int entering_variable = n;
    double var_reduced_cost = 0;

    for(size_t i = 0; i < non_basics.size(); i++)
    {   
        int var = non_basics[i];
        double reduced_cost = c(var) - y.transpose()  * A.col(var);

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

    return {entering_variable, var_reduced_cost};
    
}

std::tuple<int, bool, double> Simplex::leaving_variable(const VectorXd& d, const std::pair<int,double>& entering_variable)
{
    int leaving_variable = entering_variable.first;
    double t = std::numeric_limits<double>::infinity();
    int signal_t;
    if(entering_variable.second < -EPSILON_1)
    {
        signal_t = 1;
    }
    else
    {
        signal_t = -1;
    }

    for(int i = 0; i < m; i++)
    {
        int var = basics[i];
        double bound = -1;
        if(abs(d(i)) <= EPSILON_1) // d == 0
        {
            continue;
        }

        if(signal_t * d(i) < -EPSILON_1) // incremento negativo
        {
            bound = (x(var) - lb(var))/abs(d(i));
            if(x(var) - lb(var) <= epsilon)
            {
                bound = 0.0;
            }

        }
        else if(signal_t * d(i) > EPSILON_1) //incremento positivo
        {
            bound = (ub(var) - x(var))/abs(d(i));
            if(ub(var) - x(var) <= epsilon)
            {
                bound = 0.0;
            }
        }

        if(t - bound >= epsilon) //bound menor ou igual;
        {
            if(abs(t - bound) <= epsilon) //igual
            {
                if(var < leaving_variable)
                {
                    leaving_variable = var;
                    t = bound;
                }
            }
            else
            {
                t = bound;
                leaving_variable = var;
            }

        }
    }
    //testar variavel nao basica (parametro)

    if(t - (ub(entering_variable.first) - lb(entering_variable.first)) > epsilon)
    {
        t = ub(entering_variable.first) - lb(entering_variable.first);
        leaving_variable = entering_variable.first;
    }


    //sinal de t na perspectiva das variaves basicas
    //para nao basicas: -signal_t
    return {leaving_variable, leaving_variable != entering_variable.first, signal_t * t};
    
}

