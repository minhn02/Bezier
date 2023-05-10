#include <cassert>
#include "bezier.h"
#include <iostream>
#include <Eigen/Dense>
#include "matplotlibcpp.h"

using namespace Eigen;
namespace plt = matplotlibcpp;

void visualize_opt_bezier(long t1, long t_t, std::vector<VectorXd> points) {
    Bezier::Curve<long> curve(points, t_t - t1, t1);
    std::vector<double> x, y;

    double delta = (t_t - t1) / 1000;
    for (int i = 0; i < 1000; i++) {
        long t = t1 + delta * i;
        VectorXd p = curve.evaluate(t);
        x.push_back(t);
        y.push_back(p(0));
    }
    plt::plot(x, y);
    plt::show();
}

int main(int argc, char const *argv[])
{
    VectorXd P0(2); P0 << -0.39835453021027050768, 0;
    VectorXd P1(2); P1 << -0.39835453021027050768, 0;
    VectorXd P2(2); P2 << 0, 0;
    VectorXd P3(2); P3 << 0, 0;
    long t1 = 1683756053853286144;
    long t_t = 1683756055720939008;
    // visualize_opt_bezier(t1, t_t, {P0, P1, P2, P3});

    Bezier::Curve<long> curve({P0, P1, P2, P3}, t_t - t1, t1);
    std::cout << curve.evaluate(1683756055713785099).sum() << std::endl;

    return 0;
}
