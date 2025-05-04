#include "aux.h"

void print_vector_int(const std::vector<int>& v)
{
    std::cout << "VEC: ";
    for(size_t i = 0; i < v.size(); i++)
    {
        std::cout << v[i] << " ";
    }
    std::cout << std::endl;
}

inline void update_basis_info(Simplex& simplex, const VectorXd& d, const int entering_variable, 
    const int leaving_variable)
{

    simplex.c_b(simplex.basics_idx[leaving_variable]) = simplex.c(entering_variable);
    std::swap(simplex.basics[simplex.basics_idx[leaving_variable]], simplex.non_basics[simplex.non_basics_idx[entering_variable]]);
    std::swap(simplex.basics_idx[leaving_variable], simplex.basics_idx[entering_variable]);
    std::swap(simplex.non_basics_idx[leaving_variable], simplex.non_basics_idx[entering_variable]);
}
//TODO: tirar isso
inline void update_first_phase(Simplex& simplex, const VectorXd& d, const int entering_variable, 
    const int leaving_variable)
{
    std::cout << "CHAMOU" << std::endl;
    std::cout << "ENTER: " << entering_variable << std::endl;
    std::cout << "LEAVE: " << leaving_variable << std::endl;
    // std::cout << "C: " << simplex.c_b(simplex.basics_idx[leaving_variable]) << std::endl;
    simplex.c_b(simplex.basics_idx[leaving_variable]) = 0; //variavel ficou no bound
    simplex.c(leaving_variable) = 0; //TODO: Usar a funcao update_basis_info (mesma coisa)
    std::cout << "C_b: " << simplex.c_b(simplex.basics_idx[leaving_variable]) << std::endl;
    // std::cout << "C: " << simplex.c(simplex.basics_idx[leaving_variable]) << std::endl;
    // std::cout << "AQ" << std::endl;
    std::swap(simplex.basics[simplex.basics_idx[leaving_variable]], simplex.non_basics[simplex.non_basics_idx[entering_variable]]);
    std::swap(simplex.basics_idx[leaving_variable], simplex.basics_idx[entering_variable]);
    std::swap(simplex.non_basics_idx[leaving_variable], simplex.non_basics_idx[entering_variable]);
}

//TODO: ajeitar essas duas funcoes
void solve(const mpsReader& Data)
{
    Simplex solver(Data.A, Data.lb, Data.ub, Data.c);
    SystemSolver SS(solver.basics, solver.A);
    VectorXd init_x_b = SS.solve_initial(Data.b);
    VectorXd init_x_n = VectorXd::Zero(solver.n - solver.m);
    // solver.init_sol(init_x_b, init_x_n);
    VectorXd y;
    VectorXd d;
    double obj = 0;
    double t = 0;
    bool is_basic; // 1 sse a variavel que vai sair é basica
    // std::cout << "KD" << std::endl;
    while(true)
    {
        // std::cout << "C_b: " << solver.c_b.transpose() << std::endl;
        y = SS.solve_price(solver.c_b);
        // std::cout << "Duais: " << y.transpose() << std::endl;
        std::pair<int, double> enter_var = solver.entering_variable(y);
        if(enter_var.first == solver.n)// nenhum custo reduzido melhora
        {
            std::cout << "OTIMO: " << obj << std::endl;
            break;
            
        } 
                    
        d = SS.solve_direction(solver.A.col(enter_var.first));
        std::tuple<int, bool, double> leave_var_info = solver.leaving_variable(d, enter_var);
        t = get<2>(leave_var_info);
        if(abs(t) == std::numeric_limits<double>::infinity())
        {
            break;
        }
        is_basic = get<1>(leave_var_info);
        obj += -t * enter_var.second;
        std::cout << "OBJ: " << obj << std::endl;
        solver.update_sol(d, t, enter_var.first);

        if(is_basic)
        {
            SS.E.push_back({solver.basics_idx[get<0>(leave_var_info)], d});
            update_basis_info(solver, d, enter_var.first, get<0>(leave_var_info));   
        }
    }
}

void solve_first_phase(const mpsReader& Data)
{
    Simplex solver(Data.A, Data.lb, Data.ub, Data.c);
    SystemSolver SS(solver.basics, solver.A);

    VectorXd init_x_n = VectorXd::Zero(solver.non_basics.size());
    VectorXd Nx_n = VectorXd::Zero(solver.basics.size());
    double inf = std::numeric_limits<double>::infinity();
    VectorXd y;
    VectorXd d;
    double obj = 0;
    double t = 0;
    bool is_basic; // 1 sse a variavel que vai sair é basica
    int leave_var = -1;
    for(size_t i = 0; i < solver.non_basics.size(); i++)
    {
        // std::cout << "I: " << i << std::endl;
        int var = solver.non_basics[i];
        // std::cout << "LB: " << Data.lb[var] << std::endl;
        // std::cout << "UB: " << Data.ub[var] << std::endl;
        if(Data.lb[var] == -inf  && Data.ub[var] == inf)
        {
            init_x_n(i) = 0;
        }
        else if(Data.lb[var] == -inf)
        {
            init_x_n(i) = Data.ub[var];
        }
        else
        {
            init_x_n(i) = Data.lb[var];
        }
        Nx_n += init_x_n(i) * Data.A.col(var);
    }

    VectorXd init_x_b = SS.solve_initial(Data.b - Nx_n);
    std::cout << init_x_b.transpose() << std::endl;
    obj = solver.check_initial_infeasible(init_x_b, init_x_n);
    
    while(abs(obj) > EPSILON_1)
    {
        y = SS.solve_price(solver.c_b);
        std::pair<int, double> enter_var = solver.entering_variable(y);
        if(enter_var.first == solver.n)// nenhum custo reduzido melhora
        {
            std::cout << "OTIMO: " << obj << std::endl;
            break;
            
        } 
                    
        d = SS.solve_direction(solver.A.col(enter_var.first));
        std::tuple<int, bool, double> leave_var_info = solver.leaving_variable(d, enter_var);
        t = get<2>(leave_var_info);
        leave_var = get<0>(leave_var_info);
        if(abs(t) == std::numeric_limits<double>::infinity())
        {
            break;
        }
        is_basic = get<1>(leave_var_info);
        // std::cout << "OBJ(1): " << obj << std::endl;
        obj += -t * enter_var.second;
        std::cout << "VAR(1): " << solver.x(leave_var) << std::endl;

        solver.update_sol(d, t, enter_var.first);
        std::cout << "OBJ: " << obj << std::endl;

        std::cout << "VAR(2): " << solver.x(leave_var) << std::endl;


        if(is_basic)
        {
            SS.E.push_back({solver.basics_idx[leave_var], d});
            update_first_phase(solver, d, enter_var.first, leave_var);
            std::cout << solver.c(leave_var) << std::endl;
            std::cout << "TRUE OBJ: " << solver.c.transpose() * solver.x << std::endl; 
            if(solver.violate_bound[leave_var] == -1)
            {
                solver.lb(leave_var) = Data.lb(leave_var);
                std::cout << "VAR: " << solver.x(leave_var) << std::endl;
                std::cout << "LB: " << Data.lb(leave_var) << std::endl;
                std::cout << "OBJ ANTES: " << obj << std::endl;
                obj -= Data.lb(leave_var);
                std::cout << "OBJ AGORA: " << obj << std::endl;
                getchar();
                // obj -= Data.lb(leave_var);
                // std::cout << "OBJ agora(1): " << obj << std::endl;
            }   
            else if(solver.violate_bound[leave_var] == 1)
            {
                solver.ub(leave_var) = Data.ub(leave_var);
                std::cout << "VAR: " << solver.x(leave_var) << std::endl;
                std::cout << "UB: " << Data.ub(leave_var) << std::endl;
                std::cout << "OBJ ANTES: " << obj << std::endl;
                obj += Data.ub(leave_var);
                std::cout << "OBJ AGORA: " << obj << std::endl;
                getchar();
                // obj += Data.ub(leave_var);
                // std::cout << "OBJ agora(2): " << obj << std::endl;
            }
        }
        else
        {
            std::cout << solver.c(leave_var) << std::endl;
            std::cout << "TRUE OBJ: " << solver.c.transpose() * solver.x << std::endl; 

        }
        // getchar();

    }

    std::cout << "VIAVEL: " << obj << std::endl;




}