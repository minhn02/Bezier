#include "bezier.h"
#include "nlopt.hpp"
#include <vector>

using namespace nlopt;
using namespace Eigen;

double cost_func(const std::vector<double> &x, std::vector<double> &grad, void* f_data) {

}

int main(int argc, char const *argv[])
{
    opt optimization = opt(nlopt::LD_SLSQP, 2);
    optimization.set_min_objective(cost_func, nullptr);
    
    return 0;
}
