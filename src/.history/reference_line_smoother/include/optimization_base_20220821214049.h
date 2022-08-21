#ifndef OPTIMIZATION_BASE_H
#define OPTIMIZATION_BASE_H

#pragma once
#include <eigen3/Eigen/SparseCore>
#include <OsqpEigen/OsqpEigen.h>

#include <boost/log/trivial.hpp>

#include "../../log/log.h"

namespace optimization
{
namespace osqp_optimization
{
class OsqpOptimization
{
public: 
    OsqpOptimization() = default;
    virtual ~OsqpOptimization() {};

    void setNumOfVar(const int& num_var){num_of_var_ = num_var;}
    void setNumOfCons(const int& num_cons) {num_cons_ = num_cons;}
    void setHessian(Eigen::SparseMatrix<double>& hessian) {hessian_ = hessian;}
    void setGradient(Eigen::VectorXd& gradient) {gradient_ = gradient;}
    void setLinearMatrix(Eigen::SparseMatrix<double>& linear_matrix) {linear_matrix_=linear_matrix;}
    void setLowerBound(Eigen::VectorXd& lower_bound) {lower_bound_ = lower_bound;}
    void setUpperBound(Eigen::VectorXd& upper_bound) {upper_bound_ = upper_bound;}
    
    int getNumOfVar() const {return num_of_var_;} 
    Eigen::VectorXd getSolution() {return osqp_solution_;}
    /*
    @ return if solve OSQP problem successfully, final solution stored in osqp_solution
    */
    bool solve();
    
protected:
    int num_of_var_;
    int num_cons_;
private:
    // osqp related matrix and constraints
    Eigen::SparseMatrix<double> hessian_;
    Eigen::VectorXd gradient_;
    Eigen::SparseMatrix<double> linear_matrix_;
    Eigen::VectorXd lower_bound_;
    Eigen::VectorXd upper_bound_;

    // osqp solver
    OsqpEigen::Solver solver_;
    // solution
    Eigen::VectorXd osqp_solution_;
};
}// namespace osqp_optimization
    
} // namespace optimization


#endif