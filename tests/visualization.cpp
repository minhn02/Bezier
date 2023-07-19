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
    // put the above waypoints into the list
    std::vector<VectorXd> points(21, VectorXd(2));
    points[0] << 0.219919, 0.000000;
    points[1] << 0.154532, -0.070092;
    points[2] << 0.127777, 0.003925;
    points[3] << 0.077738, -0.014935;
    points[4] << 0.036375, 0.001309;
    points[5] << -0.009374, 0.001063;
    points[6] << -0.055585, -0.000482;
    points[7] << -0.101073, -0.000256;
    points[8] << -0.147083, -0.000595;
    points[9] << -0.193604, 0.001416;
    points[10] << -0.239333, -0.000162;
    points[11] << -0.281523, -0.000392;
    points[12] << -0.327118, 0.000564;
    points[13] << -0.376792, 0.000688;
    points[14] << -0.422793, -0.000721;
    points[15] << -0.463553, -0.001756;
    points[16] << -0.528784, -0.001316;
    points[17] << -0.562053, 0.008350;
    points[18] << -0.607742, -0.000149;
    points[19] << -0.649432, 0.000263;
    points[20] << -0.696839, 0.005503;


    Bezier::Spline<int64_t> spline(points, 1699671933, 1689731935545659642);
    int64_t t_t = 1689731935545659642 + 1699671933;

    std::printf("finished? %d\n", spline.isFinished(1689731936195756763));
    std::vector<double> x, y;
    double delta = (t_t - 1689731935545659642) / 1000;
    for (int i = 0; i < 1000; i++) {
        int64_t t = 1689731935545659642 + delta * i;
        VectorXd p = spline.evaluate(t);
        x.push_back(t);
        y.push_back(p(0));
    }
    plt::plot(x, y);
    plt::show();

    return 0;
}
