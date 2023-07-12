#include "bezier.h"
#include "nlopt.hpp"
#include "nlopt.h"
#include <vector>
#include <cmath>
#include <iostream>
#include "matplotlibcpp.h"

using namespace nlopt;
using namespace Eigen;
namespace plt = matplotlibcpp;

// define gaits
double period = 2*M_PI;
auto g1 = [](double t) {VectorXd val(2); val << std::sin(t), std::cos(t); return val;};
auto d_g1 = [](double t) {VectorXd val(2); val << std::cos(t), -std::sin(t); return val;};

auto g2 = [](double t) {VectorXd val(2); val << std::cos(t) + 10, std::sin(t) + 4; return val;};
auto d_g2 = [](double t) {VectorXd val(2); val << -std::sin(t), std::cos(t); return val;};

double t1 = M_PI;

double cost_func(const std::vector<double> &x, std::vector<double> &grad, void* f_data) {
    // x is vector {t_t, t_2}
    VectorXd P0 = g1(t1);
    VectorXd P1 = P0 + (1/3)*d_g1(t1);
    VectorXd P3 = g2(x[1]);
    VectorXd P2 = P3 - (1/3)*d_g2(x[1]);

    Bezier::Curve<double> curve({P0, P1, P2, P3}, x[0] - t1, t1);
    double acc = 0;
    double vel = 0;
    for (double i = t1; i < x[0]; i+=0.01) {
        vel += std::pow(curve.dEvaluate(i).sum(), 2);
        acc += std::pow(curve.dEvaluate(2, i).sum(), 2);
    }
    double avgAcc = acc / ((x[0] - t1)/0.01);
    double avgVel = vel / ((x[0] - t1)/0.01);
    std::printf("Cost Func, first derivative avg: %f, second derivative avg: %f, t_t %f \n", avgVel/200, avgAcc/500, x[0]);
    return avgVel/200 + avgAcc/500 + x[0];
}

int main(int argc, char const *argv[]) {
    // define optimization
    std::printf("t1: %f, period: %f \n", t1, period);

    opt optimization = opt(nlopt::LN_SBPLX, 2);
    optimization.set_min_objective(cost_func, nullptr);
    optimization.set_lower_bounds({t1, 0});
    optimization.set_upper_bounds({t1 + period/2, period});
    optimization.set_maxtime(5e-2);

    std::printf("Setting lower bounds as [%f, %f] \n", optimization.get_lower_bounds()[0], optimization.get_lower_bounds()[1]);
    std::printf("Setting upper bounds as [%f, %f] \n", optimization.get_upper_bounds()[0], optimization.get_upper_bounds()[1]);

    optimization.set_xtol_rel(1e-2);
    
    std::vector<double> guess(2);
    guess[0] = t1 + period/2; guess[1] = period/2;
    double obj_value;

    try {
        result res = optimization.optimize(guess, obj_value);
        std::printf("Optimized Vector {%f, %f}, with objective value %f \n", guess[0], guess[1], obj_value);
    } catch (std::exception &e) {
        std::printf("Optimization failed: %s \n", e.what());
        std::printf("Optimized Vector {%f, %f}, with objective value %f \n", guess[0], guess[1], obj_value);
    }

    VectorXd P0 = g1(t1);
    VectorXd P1 = P0 + (1/3)*d_g1(t1);
    VectorXd P3 = g2(guess[1]);
    VectorXd P2 = P3 - (1/3)*d_g2(guess[1]);

    Bezier::Curve<double> curve({P0, P1, P2, P3}, guess[0] - t1, t1);

    //plot output
    std::vector<double> xVals, g1Vals, g2Vals, splineVals;
    int num_points = 100;
    double time_delta = (guess[0] - t1)/(num_points);
    for (int i = 0; i < num_points; i++) {
        double x = t1 + time_delta*i;
        xVals.push_back(x);
        g1Vals.push_back(g1(x)(0));
        g2Vals.push_back(g2(x - guess[0] + guess[1])(0));
        splineVals.push_back(curve.evaluate(x)(0));
    }

    plt::plot(xVals, g1Vals, "r");
    plt::plot(xVals, g2Vals, "b");
    plt::plot(xVals, splineVals, "g");
    plt::show();

    return 0;
}
