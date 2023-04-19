#include "bezier.h"
#include "nlopt.hpp"
#include "nlopt.h"
#include <vector>
#include <cmath>
#include <iostream>

using namespace nlopt;
using namespace Eigen;

// define gaits
double period = 2*M_PI;
auto g1 = [](double t) {VectorXd val(2); val << std::sin(t), std::cos(t); return val;};
auto d_g1 = [](double t) {VectorXd val(2); val << std::cos(t), -std::sin(t); return val;};

auto g2 = [](double t) {VectorXd val(2); val << std::cos(t), std::sin(t); return val;};
auto d_g2 = [](double t) {VectorXd val(2); val << -std::sin(t), std::cos(t); return val;};

double t1 = 0;

double cost_func(const std::vector<double> &x, std::vector<double> &grad, void* f_data) {
    // x is vector {t_t, t_2}
    grad = {1, 0};
    return x[0];
}

double bezier_boundary_constraint(const std::vector<double> &x, std::vector<double> &grad, void* f_data) {
    // x is vector {t_t, t_2}
    // f_data is []
    VectorXd P0 = g1(t1);
    VectorXd P1 = P0 + (1/3)*d_g1(t1);
    VectorXd P3 = g2(x[1]);
    VectorXd P2 = P3 - (1/3)*d_g2(x[1]);

    Bezier::Curve curve({P0, P1, P2, P3}, x[0]);
    return (curve.evaluate(x[1] - t1) - g2(x[1])).sum();
}

int main(int argc, char const *argv[]) {
    // define optimization

    opt optimization = opt(nlopt::LD_SLSQP, 2);
    optimization.set_min_objective(cost_func, nullptr);
    optimization.set_lower_bounds({t1, 0});
    optimization.set_upper_bounds({t1+period, period});
    optimization.add_equality_constraint(bezier_boundary_constraint, nullptr, 1e-8);
    
    std::vector<double> guess(2);
    guess[0] = period + period/2; guess[1] = period/2;
    double obj_value;

    try {
        result res = optimization.optimize(guess, obj_value);
        std::printf("Optimized Vector {%f, %f}, with objective value %f \n", guess[0], guess[1], obj_value);
    } catch (std::exception &e) {
        std::printf("Optimization failed: %s \n", e.what());
    }

    return 0;
}
