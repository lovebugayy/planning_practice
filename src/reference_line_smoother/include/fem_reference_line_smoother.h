#ifndef FEM_REFERENCE_LINE_SMOOTHER_H
#define FEM_REFERENCE_LINE_SMOOTHER_H

#pragma once
#include <vector>
#include "optimization_base.h"

namespace optimization
{
namespace osqp_optimization
{
class FemReferenceLineSmoother : public OsqpOptimization
{
public:
    FemReferenceLineSmoother() {};
    FemReferenceLineSmoother(const std::vector<std::pair<double,double>>& reference_points);
    ~FemReferenceLineSmoother(){};

    bool init(const std::vector<std::pair<double,double>>& reference_points);

private:
    void generateHessian(Eigen::SparseMatrix<double>& hessian,const std::vector<std::pair<double,double>>& reference_points);
    void generateSmoothCostHessian(Eigen::SparseMatrix<double>& smooth_cost_hessian,const std::vector<std::pair<double,double>>& reference_points);
    void generateCloseCostHessian(Eigen::SparseMatrix<double>& close_cost_hessian,const std::vector<std::pair<double,double>>& reference_points);
    void generateDistributionCostHessian(Eigen::SparseMatrix<double>& distribution_cost_hessian,const std::vector<std::pair<double,double>>& reference_points);

    void generateCloseCostGradient( Eigen::VectorXd& close_cost_gradient,const std::vector<std::pair<double,double>>& reference_points);
    void generateLinearMatrix(Eigen::SparseMatrix<double>& linear_matrix,const std::vector<std::pair<double,double>>& reference_points);
    void generateLowerBound(Eigen::VectorXd& lower_bound,const std::vector<std::pair<double,double>>& reference_points);
    void generateUpperBound(Eigen::VectorXd& upper_bound,const std::vector<std::pair<double,double>>& reference_points);
private:
    int num_of_point_;

    double close_weight;
    double smooth_weight;
    double distribution_weight;

    double x_buffer_;
    double y_buffer_;
};
}// namespace osqp_optimization
} // namespace optimization

#endif