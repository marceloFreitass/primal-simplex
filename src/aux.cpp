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

inline void update_basis_info(Simplex& simplex, VectorXd& d, int entering_variable, 
    int leaving_variable)
{

    simplex.c_b(simplex.basics_idx[leaving_variable]) = simplex.c(entering_variable);
    std::swap(simplex.basics[simplex.basics_idx[leaving_variable]], simplex.non_basics[simplex.non_basics_idx[entering_variable]]);
    std::swap(simplex.basics_idx[leaving_variable], simplex.basics_idx[entering_variable]);
    std::swap(simplex.non_basics_idx[leaving_variable], simplex.non_basics_idx[entering_variable]);
}
void solve(mpsReader& Data, Simplex& simplex, SystemSolver& system_solver)
{

    VectorXd y;
    VectorXd d;
    double t = 0;
    size_t MAX_ETA_SIZE = 20;
    bool is_basic; // 1 sse a variavel que vai sair é basica
    simplex.init_c(Data.c);
    double obj = simplex.c.transpose() * simplex.x;
    std::cout << "OBJ: " << obj << std::endl;
    while(true)
    {
        y = system_solver.solve_price(simplex.c_b);
        std::pair<int, double> enter_var = simplex.entering_variable(y);
        if(enter_var.first == simplex.n)// nenhum custo reduzido melhora
        {
            std::cout << "OPTIMAL: " << simplex.c.transpose() * simplex.x << std::endl;
            break;
            
        } 
                    
        d = system_solver.solve_direction(simplex.A.col(enter_var.first));
        std::tuple<int, bool, double> leave_var_info = simplex.leaving_variable(d, enter_var);
        t = get<2>(leave_var_info);
        if(abs(t) == std::numeric_limits<double>::infinity())
        {
            std::cout << "UNBOUNDED" << std::endl;
            break;
        }
        is_basic = get<1>(leave_var_info);
        obj += -t * enter_var.second;
        simplex.update_sol(d, t, enter_var.first);

        if(is_basic)
        {
            system_solver.E.push_back({simplex.basics_idx[get<0>(leave_var_info)], d});
            update_basis_info(simplex, d, enter_var.first, get<0>(leave_var_info)); 
            if(system_solver.E.size() >= MAX_ETA_SIZE)
            {   
                system_solver.refactor(simplex.basics, simplex.A);
            }  
        }
        std::cout << "OBJ: " << simplex.c.transpose() * simplex.x << std::endl;

    }
}

void solve_first_phase(mpsReader& Data, Simplex& simplex, SystemSolver& system_solver)
{
    VectorXd init_x_n = VectorXd::Zero(simplex.non_basics.size());
    VectorXd Nx_n = VectorXd::Zero(simplex.basics.size());
    double inf = std::numeric_limits<double>::infinity();
    
    VectorXd y;
    VectorXd d;
    int infeasible_variables = 0;
    double t = 0;
    bool is_basic; // 1 sse a variavel que vai sair é basica
    int leave_var = -1;
    size_t MAX_ETA_SIZE = 20;

    for(size_t i = 0; i < simplex.non_basics.size(); i++)
    {
        int var = simplex.non_basics[i];
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

    VectorXd init_x_b = system_solver.solve_initial(Data.b - Nx_n);
    simplex.init_sol(init_x_b, init_x_n);
    infeasible_variables = simplex.check_infeasible(Data.lb, Data.ub);
    
    while(infeasible_variables > 0)
    {
        y = system_solver.solve_price(simplex.c_b);
        std::pair<int, double> enter_var = simplex.entering_variable(y);
        if(enter_var.first == simplex.n)// nenhum custo reduzido melhora
        {
            std::cout << "OPTIMAL: " << infeasible_variables << std::endl;
            break;
            
        } 
                    
        d = system_solver.solve_direction(simplex.A.col(enter_var.first));
        std::tuple<int, bool, double> leave_var_info = simplex.leaving_variable(d, enter_var);
        t = get<2>(leave_var_info);
        leave_var = get<0>(leave_var_info);

        if(abs(t) == std::numeric_limits<double>::infinity())
        {
            std::cout << "UNBOUNDED\n";
            break;
        }
        is_basic = get<1>(leave_var_info);

        simplex.update_sol(d, t, enter_var.first);


        if(is_basic)
        {
            system_solver.E.push_back({simplex.basics_idx[leave_var], d});
            update_basis_info(simplex, d, enter_var.first, leave_var);
            if(system_solver.E.size() >= MAX_ETA_SIZE)
            {   
                system_solver.refactor(simplex.basics, simplex.A);
            }
        }
        infeasible_variables = simplex.check_infeasible(Data.lb, Data.ub);
        std::cout << "# OUT OF BOUNDS: " << infeasible_variables << std::endl;


    }

    std::cout << "FEASIBLE: " << infeasible_variables << std::endl;
}