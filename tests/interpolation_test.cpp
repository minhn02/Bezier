#include <cassert>
#include "bezier.h"
#include <cmath>
#include <functional>
#include "matplotlibcpp.h"
#include <numeric>

using namespace Eigen;
namespace plt = matplotlibcpp;

void plot_sin_interpolation() {
    auto sin_func = [] (double t) { VectorXd ret(1); ret << std::sin(t); return ret; };
    double period = 2*M_PI;
    Bezier::Spline spline(sin_func, 100, period);

    std::vector<double> x(500);
    std::vector<double> y(500);
    std::iota(x.begin(), x.end(), 1);

    int index = 0;
    for (double t : x) {
        y[index] = spline.evaluate(t).sum();
        index++;
    }

    plt::plot(x, y);
    plt::show();
}

// void test_sin_interpolation() {
//     auto sin_func = [] (double t) { VectorXd ret(1); ret << std::sin(t); return ret; };
//     double period = 2*M_PI;
//     Bezier::Spline spline(sin_func, period, 100);

//     double MSE = 0;
//     double delta = 0.01;
//     for (double i = 0; i < period; i+=delta) {
//         MSE += std::pow(spline.evaluate(i).sum() - std::sin(i), 2);
//     }
//     MSE = MSE / (period/delta);
//     assert(MSE <= 0.01);
// }

// void test_sin_interpolation_derivative() {
//     auto sin_func = [] (double t) { VectorXd ret(1); ret << std::sin(t); return ret; };
//     double period = 2*M_PI;
//     Bezier::Spline spline(sin_func, period, 100);

//     double MSE = 0;
//     double delta = 0.01;
//     for (double i = 0; i < period; i+=delta) {
//         MSE += std::pow(spline.dEvaluate(i).sum() - std::cos(i), 2);
//     }
//     MSE = MSE / (period/delta);
//     assert(MSE <= 0.01);
// }

// void test_abs_interpolation() {
//     auto abs_func = [] (double t) { VectorXd ret(1); ret << std::abs(t); return ret; };
//     double period = 8;
//     //TODO add where spline starts
//     Bezier::Spline spline(abs_func, period, -period/2, 100);

//     double delta = 0.01;
//     double MSE = 0;
//     for (double i = -period/2; i < period/2; i+=delta) {
//         MSE += std::pow(spline.evaluate(i).sum() - std::abs(i), 2);
//     }
//     MSE = MSE / (period/delta);
//     assert(MSE <= 0.01);
// }

// void test_sin_multidim_interpolation() {
//     auto sin_func = [] (double t) { VectorXd ret(2); ret << std::sin(t), std::cos(t); return ret; };
//     double period = 2*M_PI;
//     Bezier::Spline spline(sin_func, period, 100);

//     double MSE = 0;
//     double delta = 0.01;
//     for (double i = 0; i < period; i+=delta) {
//         VectorXd comp(2); comp << std::sin(i), std::cos(i);
//         MSE += std::pow(spline.evaluate(i).sum() - comp.sum(), 2);
//     }
//     MSE = MSE / (period/delta);
//     assert(MSE <= 0.01);
// }

int main(int argc, char const *argv[]) {
    plot_sin_interpolation();
    // test_sin_interpolation();
    // test_sin_interpolation_derivative();
    // test_abs_interpolation();
    // test_sin_multidim_interpolation();
}
