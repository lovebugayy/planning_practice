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

    std::vector<double> x,y;
    for(auto i=0;i<5;i++)
    {
        x.emplace_back(solution(i));
        y.emplace_back(solution(i+5));
        std::cout << "("<<solution(i)<<","<<solution(i+5)<<")"<< std::ends;
    }

    // plt::plot(x,y,"r-");
    
    std::vector<double> t(1000);
    std::vector<double> x(t.size());

    for(size_t i = 0; i < t.size(); i++) {
        t[i] = i / 100.0;
        x[i] = sin(2.0 * M_PI * 1.0 * t[i]);
    }

    plt::xkcd();
    plt::plot(t, x);
    plt::title("AN ORDINARY SIN WAVE");
    plt::save("xkcd.png");
    return 0;
}
