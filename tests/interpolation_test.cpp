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
    Bezier::Spline spline(sin_func, 500, period, -period/2);

    double n = 3000;
    double delta = period/n;
    std::vector<double> x(n);
    std::vector<double> y(n);

    int index = 0;
    for (double t = -period/2; t <= period/2; t+=delta) {
        x[index] = t;
        y[index] = spline.evaluate(t).sum();
        index++;
    }

    plt::plot(x, y);
    plt::show();
}

void test_sin_interpolation() {
    auto sin_func = [] (double t) { VectorXd ret(1); ret << std::sin(t); return ret; };
    double period = 2*M_PI;
    Bezier::Spline spline(sin_func, 500, period);

    double MSE = 0;
    double delta = 0.01;
    int n = period/delta;

    std::vector<double> x(n);
    std::vector<double> ySpline(n);
    std::vector<double> ySin(n);
    std::vector<double> MSEs(n);

    int index = 0;
    for (double t = 0; t <= period; t+=delta) {
        x[index] = t;
        ySpline[index] = spline.evaluate(t).sum();
        ySin[index] = std::sin(t);
        MSE += std::pow(ySpline[index] - ySin[index], 2);
        MSEs[index] = MSE/index;
        index++;
    }
    MSE = MSE / (double) n;

    // plt::plot(x, ySin, "r");
    // plt::plot(x, ySpline, "g");
    // plt::plot(x, MSEs, "b");
    // plt::show();
    assert(MSE <= 0.01);
}

void test_sin_interpolation_derivative() {
    auto sin_func = [] (double t) { VectorXd ret(1); ret << std::sin(t); return ret; };
    double period = 2*M_PI;
    Bezier::Spline spline(sin_func, 500, period);

    double MSE = 0;
    double delta = .01;
    int n = period/delta;
    std::vector<double> x(n);
    std::vector<double> dCurve(n);
    std::vector<double> curve(n);

    int index = 0;
    for (double t = 0; t <= period; t+=delta) {
        MSE += std::pow(spline.dEvaluate(t).sum() - std::cos(t), 2);

        x[index] = t;
        dCurve[index] = spline.dEvaluate(t).sum();
        curve[index] = spline.evaluate(t).sum();
        index++;
    }
    MSE = MSE / (double)n;

    plt::plot(x, dCurve, "b");
    plt::plot(x, curve, "r");
    plt::show();

    assert(MSE <= 0.01);
}

void test_abs_interpolation() {
    auto abs_func = [] (double t) { VectorXd ret(1); ret << std::abs(t); return ret; };
    double period = 10;
    Bezier::Spline spline(abs_func, 500, period, -period/2);

    double MSE = 0;
    double delta = 0.01;
    double n = period/delta;

    for (double t = -period/2; t <= period/2; t+=delta) {
        MSE += std::pow(spline.evaluate(t).sum() - std::abs(t), 2);
    }
    MSE = MSE / n;
    assert(MSE <= 0.01);
}

void test_sin_multidim_interpolation() {
    auto sin_func = [] (double t) { VectorXd ret(2); ret << std::sin(t), std::cos(t); return ret; };
    double period = 2*M_PI;
    Bezier::Spline spline(sin_func, 500, period);

    double MSE = 0;
    double delta = 0.01;
    for (double i = 0; i <= period; i+=delta) {
        VectorXd comp(2); comp << std::sin(i), std::cos(i);
        MSE += std::pow(spline.evaluate(i).sum() - comp.sum(), 2);
    }
    MSE = MSE / (period/delta);
    assert(MSE <= 0.01);
}

int main(int argc, char const *argv[]) {
    // plot_linear_interpolation();
    test_sin_interpolation();
    test_sin_interpolation_derivative();
    // test_abs_interpolation();
    // test_sin_multidim_interpolation();
}
