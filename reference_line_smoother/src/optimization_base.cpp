#include "optimization_base.h"

namespace optimization
{
namespace osqp_optimization
{
bool OsqpOptimization::solve()
{
    //settings 
    std::cout << "setting for solver base" << std::endl;
    // solver_.settings()->setVerbosity(false);
    solver_.settings()->setWarmStart(true);

    //set initial data of OSQP sovler
    std::cout << "num_of_var_: "<<num_of_var_ <<", num_cons_"<<num_cons_;
    std::cout << "init for solver base" << std::endl;
    
    solver_.data()->setNumberOfVariables(num_of_var_);
    std::cout << "init nums constraints" << std::endl;
    solver_.data()->setNumberOfConstraints(num_cons_);

    std::cout << "matrix and constraints for solver base" << std::endl;
    if(!solver_.data()->setHessianMatrix(hessian_)) return false;
    if(!solver_.data()->setGradient(gradient_)) return false;
    if(!solver_.data()->setLinearConstraintsMatrix(linear_matrix_)) return false;
    if(!solver_.data()->setLowerBound(lower_bound_)) return false;
    if(!solver_.data()->setUpperBound(upper_bound_)) return false;

    // initial the solver
    if(!solver_.initSolver()) return false;

    // solve osqp problem
    if(!solver_.solve()) return false;

    // get solution
    osqp_solution_ = solver_.getSolution();
    return true;
}
}// namespace osqp_optimization
    
} // namespace optimization