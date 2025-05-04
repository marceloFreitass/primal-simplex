#include "Simplex.h"

Simplex::Simplex(const MatrixXd& A, const VectorXd& lb, const VectorXd& ub, const VectorXd& c)
{

    this->A = A.sparseView();
    this->lb = lb;
    this->ub = ub;
    n = A.cols(); //ja ta com de folga, numero real de variaives = n - m
    m = A.rows();

    // this->c = c; TODO
    this->c = VectorXd::Zero(n);

    basics.reserve(m);
    non_basics.reserve(n - m);
    basics_idx = std::vector<int>(n, -1);
    non_basics_idx = std::vector<int>(n, -1);
    violate_bound = std::vector<int>(n, 0);
    c_b = VectorXd::Zero(m);
    // std::cout << c_b << std::endl;
    // std::cout << "N: " << n << std::endl;
    // std::cout << "M: " << m << std::endl;
    x = VectorXd::Zero(n);
    // x_b = VectorXd::Zero(m); // TODO: CONSEGUIR SOLUCAO INICIAL DE FORMA ALGORITMICA
    // x_n = VectorXd::Zero(n - m);
    //forma B = I;
    for(int i = n - m; i < n; i++)
    {
        // std::cout << "ind: " << i - (n - m) << "\n";
        basics.push_back(i);
        c_b(i - (n - m)) = c(basics[i - (n - m)]);
        basics_idx[i] = (i - (n - m));
        // aux++;
        // std::cout << c_b << std::endl;
    }
    for(int i = 0; i < n - m; i++)
    {
        non_basics.push_back(i);
        non_basics_idx[i] = i;
        // this->c[non_basics[i]] = 0;
    }

    // print_vector_int(basics);
    // print_vector_int(non_basics);

    // std::cout << c_b << std::endl;
}

double Simplex::check_initial_infeasible(const VectorXd& init_x_b, const VectorXd& init_x_n)
{
    
    for(int i = 0; i < m; i++)
    {
        x(basics[i]) = init_x_b(i);
    }
    for(int i = 0; i < n - m; i++)
    {
        x(non_basics[i]) = init_x_n(i);
    }
    //FAZENDO COM QUE C_b seja -1 ou 1 e C_n = 0 (inicializado)
    //TODO: REGULARIZAR DEPOIS
    double inf = std::numeric_limits<double>::infinity();
    double total_inf = 0;
    for(int i = 0; i < m; i++)
    {
    
        if(lb(basics[i]) - x(basics[i]) > EPSILON_1)
        {
            std::cout << "VAR: " << basics[i] << std::endl;
            std::cout << "VIO L: " << lb(basics[i]) - x(basics[i]) << std::endl;
            std::cout << "LB: " <<  lb(basics[i]) << std::endl;
            std::cout << "X: " << x(basics[i]) << std::endl;
            c_b(i) = 1;
            c(basics[i]) = 1;
            total_inf += x(basics[i]);
            lb(basics[i]) = -inf;
            violate_bound[basics[i]] = -1;
        }
        else if(x(basics[i]) - ub(basics[i]) > EPSILON_1)
        {
            std::cout << "VAR: " << basics[i] << std::endl;
            std::cout << "VIO U: " << x(basics[i]) - ub(basics[i]) << std::endl;
            std::cout << "UB: " <<  ub(basics[i]) << std::endl;
            std::cout << "X: " << x(basics[i]) << std::endl;
            c_b(i) = -1;
            c(basics[i]) = -1;
            total_inf -= x(basics[i]);
            ub(basics[i]) = inf;
            violate_bound[basics[i]] = 1;
        }
        else
        {
            c_b(i) = 0;
        }
    }
    
    
    std::cout << "INF: " << total_inf << std::endl;
    return total_inf;



    
}

std::pair<int,double> Simplex::entering_variable(const VectorXd& y)
{
    int entering_variable = n;
    std::cout << "N: " << entering_variable << std::endl;
    double var_reduced_cost = 0;
    std::cout << "Duals: \n" << y.transpose() << std::endl;
    for(size_t i = 0; i < non_basics.size(); i++)
    {   
        int var = non_basics[i];
        double reduced_cost = c(var) - y.transpose()  * A.col(var);
        double test = 0;
        VectorXd col = A.col(var);
        // std::cout << "Col: " << col.transpose() << std::endl;
        // for(int j = 0; j < y.size(); j++)
        // {
        //     test += y(j) * col(j);
        //     if(y(j) * col(j) > EPSILON_1)
        //     {
        //         std::cout << "J: " << j << std::endl;
        //         std::cout << "Teste: " << test << std::endl;
        //         getchar();

        //     }
        // }

        // std::cout << "RC: " << reduced_cost << std::endl;
        // getchar();
        // std::cout << "COLUNA: \n";
        // std::cout << A.col(var) << std::endl;
        
        if(reduced_cost > EPSILON_1)
        {
            // std::cout << "A!" << std::endl;
            // getchar();
        }
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

std::tuple<int, bool, double> Simplex::leaving_variable(const VectorXd& d, const std::pair<int,double>& entering_variable)
{
    bool basic = false;
    int leaving_variable = entering_variable.first;
    // std::cout << "ENTROU: " << leaving_variable << std::endl;
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
    // std::cout << "sinal t: " << signal_t << std::endl;
    // std::cout << "A!" << std::endl;
    for(int i = 0; i < m; i++)
    {
        int var = basics[i];
        // if(var == 6593)
        // {
        //     std::cout << "AQQQQQQQQQ" << std::endl;
        //     std::cout << "INCREMENTO: " << signal_t * t * d(i) << std::endl;
        //     std::cout << "T: " << t << std::endl;
        //     std::cout << "D: " << d(i) << std::endl;
        // }
        double bound = -1; //checar se pode dar ruim
        if(abs(d(i)) <= EPSILON_1) // d == 0
        {
            continue;
        }
        if(signal_t * d(i) < -EPSILON_1) // incremento negativo
        {
            // std::cout << "CASO 1" << std::endl;s
            bound = (x(var) - lb(var))/abs(d(i));

        }
        else if(signal_t * d(i) > EPSILON_1) //incremento positivo
        {
            // std::cout << "CASO 2: " << std::endl;
            bound = (ub(var) - x(var))/abs(d(i));
        }
        if(t - bound > EPSILON_1) //bound menor;
        {
            t = bound;
            leaving_variable = var;
            basic = true;
        }
        else if(abs(t - bound) <= EPSILON_1) //bound igual
        {
            if(var < leaving_variable)
            {
                leaving_variable = var;
            }
        }
    }
    //testar variavel nao basica (parametro)


    if(!basic) //a variavel nao basica so trouca de bound
    {
        if(signal_t == 1) //xj(t) = xj - t
        {
            t = x(leaving_variable) - lb(leaving_variable);
        }
        else //xj(t) = xj + t
        {
            t = ub(leaving_variable) - x(leaving_variable);
        }    
    }

    std::cout << "T: " << t << std::endl;
    std::cout << "VAI SAIR: " << leaving_variable << std::endl;
    std::cout << "BASICA: " << basic << std::endl;

    //sinal de t na perspectiva das variaves basicas
    //para nao basicas: -signal_t
    return {leaving_variable, basic, signal_t * t};
    
}

