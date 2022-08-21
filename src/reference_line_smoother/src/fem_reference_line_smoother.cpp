#include "fem_reference_line_smoother.h"

namespace optimization
{
namespace osqp_optimization
{
FemReferenceLineSmoother::FemReferenceLineSmoother(const std::vector<std::pair<double,double>>& reference_points)
    :OsqpOptimization()
{
    std::cout << "set related parameter." << std::endl; 
    num_of_point_ = reference_points.size();
    num_of_var_ = 2*num_of_point_;

    close_weight = 1;
    smooth_weight = 1;
    distribution_weight = 1;

    x_buffer_ = 0.2;
    y_buffer_ = 0.2;
    
     
    if(!init(reference_points));
    {
        std::cout << "input reference points is less than 3 points." << std::endl;
        return ;
    }
    std::cout << "init smoother." << std::endl;
}

bool FemReferenceLineSmoother::init(const std::vector<std::pair<double,double>>& reference_points)
{
    if(reference_points.size()<3)
    {
        return false;
    }
    Eigen::SparseMatrix<double> hessian,linear_matrix;
    Eigen::VectorXd gradient,lower_bound,upper_bound;

    int nums_cons = reference_points.size() * 2;
    std::cout << "calculate matrix and constraints." << std::endl; 
    generateHessian(hessian,reference_points);
    std::cout << "hessian size: " << hessian.outerSize() << std::endl;
    for(auto i=0;i<hessian.outerSize();i++)
    {
        for(Eigen::SparseMatrix<double>::InnerIterator it(hessian,i);it;++it)
        {
            std::cout << "(" << it.row() <<"," << it.col() <<") = "<<it.value() << ",";
        }
        std::cout<<std::endl;
    }
    generateCloseCostGradient(gradient,reference_points);
    std::cout<<"gradients: "<<std::endl;
    for(auto i=0;i<gradient.size();i++)
    {
        std::cout << gradient(i);
        
        std::cout<<std::endl;
    }
    generateLinearMatrix(linear_matrix,reference_points);
    std::cout << "linear_matrix size: " << hessian.outerSize() << std::endl;
    for(auto i=0;i<linear_matrix.outerSize();i++)
    {
        for(Eigen::SparseMatrix<double>::InnerIterator it(linear_matrix,i);it;++it)
        {
            std::cout << "(" << it.row() <<"," << it.col() <<") = "<<it.value() << ",";
        }
        std::cout<<std::endl;
    }
    generateLowerBound(lower_bound,reference_points);
     std::cout<<"lower_bound: "<<std::endl;
    for(auto i=0;i<lower_bound.size();i++)
    {
        std::cout << lower_bound(i);
        
        std::cout<<std::endl;
    }
    generateUpperBound(upper_bound,reference_points);
    std::cout<<"upper_bound: "<<std::endl;
    for(auto i=0;i<upper_bound.size();i++)
    {
        std::cout << upper_bound(i);
        
        std::cout<<std::endl;
    }

    std::cout << "set matrix and constraints." << std::endl; 
    setNumOfCons(nums_cons);
    setHessian(hessian);
    setGradient(gradient);
    setLinearMatrix(linear_matrix);
    setLowerBound(lower_bound);
    setUpperBound(upper_bound);
    std::cout << "init success." << std::endl; 
    return true;
}

void FemReferenceLineSmoother::generateHessian(Eigen::SparseMatrix<double>& hessian,const std::vector<std::pair<double,double>>& reference_points)
{
    Eigen::SparseMatrix<double> smooth_cost_hessian,close_cost_hessian,distribution_cost_hessian;
    hessian.resize(num_of_var_,num_of_var_);
    std::cout << "resize success." << std::endl;
    generateSmoothCostHessian(smooth_cost_hessian,reference_points);
    for(auto i=0;i<smooth_cost_hessian.outerSize();i++)
    {
        for(Eigen::SparseMatrix<double>::InnerIterator it(smooth_cost_hessian,i);it;++it)
        {
            std::cout << "(" << it.row() <<"," << it.col() <<") = "<<it.value() << ",";
        }
        std::cout<<std::endl;
    }
    std::cout << "generateSmoothCostHessian success." << std::endl;
    generateCloseCostHessian(close_cost_hessian,reference_points);
    for(auto i=0;i<close_cost_hessian.outerSize();i++)
    {
        for(Eigen::SparseMatrix<double>::InnerIterator it(close_cost_hessian,i);it;++it)
        {
            std::cout << "(" << it.row() <<"," << it.col() <<") = "<<it.value() << ",";
        }
        std::cout<<std::endl;
    }
    std::cout << "generateCloseCostHessian success." << std::endl;
    generateDistributionCostHessian(distribution_cost_hessian,reference_points);
    for(auto i=0;i<distribution_cost_hessian.outerSize();i++)
    {
        for(Eigen::SparseMatrix<double>::InnerIterator it(distribution_cost_hessian,i);it;++it)
        {
            std::cout << "(" << it.row() <<"," << it.col() <<") = "<<it.value() << ",";
        }
        std::cout<<std::endl;
    }
    std::cout << "generateDistributionCostHessian success." << std::endl;
    // Notice: hessian should be 2 times of sum, because of the form of the osqp problem form
    hessian = 2*(smooth_weight*smooth_cost_hessian + close_weight*close_cost_hessian + distribution_weight*distribution_cost_hessian);
    return ;
}

void FemReferenceLineSmoother::generateSmoothCostHessian(Eigen::SparseMatrix<double>& smooth_cost_hessian,const std::vector<std::pair<double,double>>& reference_points)
{
    smooth_cost_hessian.resize(num_of_var_,num_of_var_);
    std::cout << "resize success: " <<  num_of_var_ <<std::endl;
    std::cout << "num of point: " <<  num_of_point_ <<std::endl;
    for(auto i=0;i<num_of_point_;i++)
    {
        //for the first row of x and y
        if(i==0)
        {
            std::cout << "first row";
            smooth_cost_hessian.insert(i,0) = 1.0;
            smooth_cost_hessian.insert(i,1) = -2.0;
            smooth_cost_hessian.insert(i,2) = 1.0;
            smooth_cost_hessian.insert(num_of_point_+i,num_of_point_) = 1.0;
            smooth_cost_hessian.insert(num_of_point_+i,num_of_point_+1) = -2.0;
            smooth_cost_hessian.insert(num_of_point_+i,num_of_point_+2) = 1.0;
        }else if(i==1)
        {
            if(num_of_point_ == 3)
            {
                smooth_cost_hessian.insert(i,0) = -2.0;
                smooth_cost_hessian.insert(i,1) = 4.0;
                smooth_cost_hessian.insert(i,2) = -2.0;
                smooth_cost_hessian.insert(num_of_point_+i,num_of_point_) = -2.0;
                smooth_cost_hessian.insert(num_of_point_+i,num_of_point_+1) = 4.0;
                smooth_cost_hessian.insert(num_of_point_+i,num_of_point_+2) = -2.0;
            }else{
                smooth_cost_hessian.insert(i,0) = -2.0;
                smooth_cost_hessian.insert(i,1) = 5.0;
                smooth_cost_hessian.insert(i,2) = -4.0;
                smooth_cost_hessian.insert(i,3) = 1.0;
                smooth_cost_hessian.insert(num_of_point_+i,num_of_point_) = -2.0;
                smooth_cost_hessian.insert(num_of_point_+i,num_of_point_+1) = 5.0;
                smooth_cost_hessian.insert(num_of_point_+i,num_of_point_+2) = -4.0;
                smooth_cost_hessian.insert(num_of_point_+i,num_of_point_+3) = 1.0;
            }
            
        }else if(i==num_of_point_-1)
        {
            std::cout << "last two row";
            smooth_cost_hessian.insert(i,num_of_point_-1) = 1.0;
            smooth_cost_hessian.insert(i,num_of_point_-2) = -2.0;
            smooth_cost_hessian.insert(i,num_of_point_-3) = 1.0;
            smooth_cost_hessian.insert(num_of_point_+i,num_of_var_-1) = 1.0;
            smooth_cost_hessian.insert(num_of_point_+i,num_of_var_-2) = -2.0;
            smooth_cost_hessian.insert(num_of_point_+i,num_of_var_-3) = 1.0;
        }else if(i==num_of_point_-2)
        {
            std::cout << "last one  row";
            smooth_cost_hessian.insert(i,num_of_point_-1) = -2.0;
            smooth_cost_hessian.insert(i,num_of_point_-2) = 5.0;
            smooth_cost_hessian.insert(i,num_of_point_-3) = -4.0;
            smooth_cost_hessian.insert(i,num_of_point_-4) = 1.0;
            smooth_cost_hessian.insert(num_of_point_+i,num_of_var_-1) = -2.0;
            smooth_cost_hessian.insert(num_of_point_+i,num_of_var_-2) = 5.0;
            smooth_cost_hessian.insert(num_of_point_+i,num_of_var_-3) = -4.0;
            smooth_cost_hessian.insert(num_of_point_+i,num_of_var_-4) = 1.0;
        }else{
            std::cout << "else row";
            smooth_cost_hessian.insert(i,i-2) = 1.0;
            smooth_cost_hessian.insert(i,i-1) = -4.0;
            smooth_cost_hessian.insert(i,i) = 6.0;
            smooth_cost_hessian.insert(i,i+1) = -4.0;
            smooth_cost_hessian.insert(i,i+2) = 1.0;
            smooth_cost_hessian.insert(num_of_point_+i,num_of_point_+i-2) = 1.0;
            smooth_cost_hessian.insert(num_of_point_+i,num_of_point_+i-1) = -4.0;
            smooth_cost_hessian.insert(num_of_point_+i,num_of_point_+i) = 6.0;
            smooth_cost_hessian.insert(num_of_point_+i,num_of_point_+i+1) = -4.0;
            smooth_cost_hessian.insert(num_of_point_+i,num_of_point_+i+2) = 1.0;
        }
    }
    return ;
}

void FemReferenceLineSmoother::generateCloseCostHessian(Eigen::SparseMatrix<double>& close_cost_hessian,const std::vector<std::pair<double,double>>& reference_points)
{
    close_cost_hessian.resize(num_of_var_,num_of_var_);
    for(auto i=0;i<num_of_var_;i++)
    {
        close_cost_hessian.insert(i,i) = 1.0;
    }
    return ;
}

void FemReferenceLineSmoother::generateDistributionCostHessian(Eigen::SparseMatrix<double>& distribution_cost_hessian,const std::vector<std::pair<double,double>>& reference_points)
{
    distribution_cost_hessian.resize(num_of_var_,num_of_var_);
    for(auto i=0;i<num_of_point_;i++)
    {
        //for the first row of x and y
        if(i==0)
        {
            distribution_cost_hessian.insert(i,0) = 1.0;
            distribution_cost_hessian.insert(i,1) = -1.0;
            distribution_cost_hessian.insert(num_of_point_+i,num_of_point_) = 1.0;
            distribution_cost_hessian.insert(num_of_point_+i,num_of_point_+1) = -1.0;
        }else if(i == num_of_point_-1){
            distribution_cost_hessian.insert(i,i) = 1.0;
            distribution_cost_hessian.insert(i,i-1) = -1.0;
            distribution_cost_hessian.insert(num_of_point_+i,num_of_point_+i) = 1.0;
            distribution_cost_hessian.insert(num_of_point_+i,num_of_point_+i-1) = -1.0;
        }else{
            distribution_cost_hessian.insert(i,i-1) = -1.0;
            distribution_cost_hessian.insert(i,i) = 2.0;
            distribution_cost_hessian.insert(i,i+1) = -1.0;
            distribution_cost_hessian.insert(num_of_point_+i,num_of_point_+i-1) = -1.0;
            distribution_cost_hessian.insert(num_of_point_+i,num_of_point_+i) = 2.0;
            distribution_cost_hessian.insert(num_of_point_+i,num_of_point_+i+1) = -1.0;
        }
    }
    return ;
}

void FemReferenceLineSmoother::generateCloseCostGradient(Eigen::VectorXd& close_cost_gradient,const std::vector<std::pair<double,double>>& reference_points)
{
    close_cost_gradient.resize(num_of_var_);
    for(auto i=0;i<num_of_point_;i++)
    {
        close_cost_gradient(i) = -2*reference_points[i].first;
        close_cost_gradient(i+num_of_point_) = -2*reference_points[i].second;
    }
    return ;
}

void FemReferenceLineSmoother::generateLinearMatrix(Eigen::SparseMatrix<double>& linear_matrix,const std::vector<std::pair<double,double>>& reference_points)
{
    linear_matrix.resize(num_of_var_,num_of_var_);
    for(auto i=0;i<num_of_var_;i++)
    {
        linear_matrix.insert(i,i) = 1.0;
    }
    return ;
}

void FemReferenceLineSmoother::generateLowerBound(Eigen::VectorXd& lower_bound,const std::vector<std::pair<double,double>>& reference_points)
{
    lower_bound.resize(num_of_var_);
    for(auto i=0;i<num_of_point_;i++)
    {
        // notice: first point and last point should be hard constraints
        if(i == 0 || i==num_of_point_-1)
        {
            lower_bound(i) = reference_points[i].first;
            lower_bound(i+num_of_point_) = reference_points[i].second;
        }else{
            lower_bound(i) = reference_points[i].first - x_buffer_;
            lower_bound(i+num_of_point_) = reference_points[i].second - y_buffer_;
        }
        
    }
}

void FemReferenceLineSmoother::generateUpperBound(Eigen::VectorXd& upper_bound,const std::vector<std::pair<double,double>>& reference_points)
{
    // notice: first point and last point should be hard constraints
    upper_bound.resize(num_of_var_);
    for(auto i=0;i<num_of_point_;i++)
    {
        if(i == 0 || i==num_of_point_-1)
        {
            upper_bound(i) = reference_points[i].first;
            upper_bound(i+num_of_point_) = reference_points[i].second;
        }else{
            upper_bound(i) = reference_points[i].first + x_buffer_;
            upper_bound(i+num_of_point_) = reference_points[i].second + y_buffer_;
        }
        
    }
}

}// namespace osqp_optimization
    
} // namespace optimization
