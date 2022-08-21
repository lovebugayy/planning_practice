#include "fem_reference_line_smoother.h"
#include

using namespace optimization::osqp_optimization;

int main(int argc,char* argv[])
{
    
    std::vector<std::pair<double, double>> ref{{1.0,2.5},{2.2,3.0},{3.3,4.1},{4.4,5.2},{5.2,6.0}};
    // std::cout << "create ref points" << std::endl;
    std::shared_ptr<FemReferenceLineSmoother> smooth_ptr = std::make_shared<FemReferenceLineSmoother>(ref);
    if(!smooth_ptr->solve())
    {
        std::cout << "smooth failed" << std::endl;
        return 1;
    }
    std::cout << "smooth success" << std::endl;
    Eigen::VectorXd solution = smooth_ptr->getSolution();
    for(auto i=0;i<5;i++)
    {
        std::cout << "("<<solution(i)<<","<<solution(i+5)<<")"<< std::ends;
    }
    return 0;
}
