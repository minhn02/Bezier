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
    int sample_n = 20;

    Bezier::Spline spline(sin_func, sample_n, period);

    std::vector<double> t_samples(sample_n);
    std::vector<double> y_samples(sample_n);

    double sample_delta = period/(double)(sample_n-1);
    double t_sample = 0;
    for (int i = 0; i < sample_n; i++) {
        t_samples[i] = t_sample;
        y_samples[i] = spline.evaluate(t_sample).sum();

        t_sample += sample_delta;
    }

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
        MSE += std::pow(ySpline[index] - ySin[index], 2) / (double) n;
        MSEs[index] = MSE;
        index++;
    }

    // plt::plot(x, ySin, "r");
    // plt::plot(x, ySpline, "g");
    // plt::plot(x, MSEs, "b");
    // plt::scatter(t_samples, y_samples, 10);
    // plt::show();
    assert(MSE <= 0.01);
}

void test_sin_interpolation_derivative() {
    auto sin_func = [] (double t) { VectorXd ret(1); ret << std::sin(t); return ret; };
    double period = 6;
    int sample_n = 7;
    Bezier::Spline spline(sin_func, sample_n, period);

    double MSE = 0;
    int n = 1000;
    double delta = period/(double)(n-1);
    std::vector<double> x(n);
    std::vector<double> dCurve(n);
    std::vector<double> curve(n);
    std::vector<double> dFunc(n);
    std::vector<double> t_samples(sample_n);
    std::vector<double> y_samples(sample_n);

    double t_sample = 0;
    for (int i = 0; i < sample_n; i++) {
        t_samples[i] = t_sample;
        y_samples[i] = sin_func(t_sample).sum();

        t_sample += period/(double)(sample_n-1);
    }

    int index = 0;
    for (double t = 0; t <= period; t+=delta) {
        MSE += std::pow(spline.dEvaluate(t).sum() - std::cos(t), 2);

        x[index] = t;
        dCurve[index] = spline.dEvaluate(t).sum();
        curve[index] = spline.evaluate(t).sum();
        dFunc[index] = std::cos(t);
        index++;
    }
    MSE = MSE / (double)n;

    plt::named_plot("dCurve", x, dCurve, "b");
    plt::named_plot("curve", x, curve, "r");
    plt::named_plot("dFunc", x, dFunc, "g");
    plt::scatter(t_samples, y_samples, 10);
    plt::legend();
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
    test_abs_interpolation();
    test_sin_multidim_interpolation();
}
