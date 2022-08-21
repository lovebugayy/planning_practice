#include <vector>

#include "fem_reference_line_smoother.h"
#include "matplotlibcpp.h"

using namespace optimization::osqp_optimization;
namespace plt = matplotlibcpp;

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

    std::vector<double> x1,y;
    for(auto i=0;i<5;i++)
    {
        x1.emplace_back(solution(i));
        y.emplace_back(solution(i+5));
        std::cout << "("<<solution(i)<<","<<solution(i+5)<<")"<< std::ends;
    }

    plt::plot(x,y);
    plt::show();
    
   
    return 0;
}
