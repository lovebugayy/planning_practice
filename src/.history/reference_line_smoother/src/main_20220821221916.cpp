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
    std::vector<std::vector<double>> x, y, z;
    for (double i = -5; i <= 5;  i += 0.25) {
        std::vector<double> x_row, y_row, z_row;
        for (double j = -5; j <= 5; j += 0.25) {
            x_row.push_back(i);
            y_row.push_back(j);
            z_row.push_back(::std::sin(::std::hypot(i, j)));
        }
        x.push_back(x_row);
        y.push_back(y_row);
        z.push_back(z_row);
    }

    plt::plot_surface(x, y, z);
    plt::show();
    return 0;
}
