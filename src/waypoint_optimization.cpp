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
int n_timesteps = 20;
double velocity_constraint = 0.67;

struct VelocityConstraintData {
    int i;
    VectorXd x0;
    int num_waypoints;
    double velocity_bound;
};

double cost_func(const std::vector<double> &x, std::vector<double> &grad, void* f_data) {
    double obj_value = 0;
    Eigen::VectorXd x0(2); x0 = g1(t1);
    Eigen::VectorXd xf(2); xf = g2(x[0]);

    // minimize distance of points from boundary positions
    Eigen::VectorXd final_waypoint(2); final_waypoint << x[2*n_timesteps], x[2*n_timesteps+1];
    double final_diff = (xf - final_waypoint).norm();
    obj_value += final_diff;

    // penalize points from linear trajectory
    Eigen::VectorXd position_delta = (xf - x0)/n_timesteps;
    for (int i = 1; i < n_timesteps+1; i++) {
        Eigen::VectorXd waypoint(2); waypoint << x[2*i], x[2*i+1];
        Eigen::VectorXd linear_waypoint = x0 + position_delta*i;
        obj_value += (linear_waypoint - waypoint).norm();
    }

    // penalize moving steering joint
    for (int i = 1; i < n_timesteps; i++) {
        obj_value += std::abs(x[2*(i+1)] - x[2*i]);
    }

    return obj_value;
}


double steering_velocity_constraint(unsigned n, const double *x, double *grad, void *data) {
    VelocityConstraintData *d = (VelocityConstraintData*)data;
    int i = d->i;
    VectorXd x0 = d->x0;
    int num_waypoints = d->num_waypoints;
    double velocity_bound = d->velocity_bound;
    if (i == 0) {
        return std::abs((x[2] - x0(0))/(x[1]/num_waypoints)) - velocity_bound;
    } else {
        return std::abs((x[(i+1)*2] - x[2*i])/(x[1]/num_waypoints)) - velocity_bound;
    }
}
double bogie_velocity_constraint(unsigned n, const double *x, double *grad, void *data) {
    VelocityConstraintData *d = (VelocityConstraintData*)data;
    int i = d->i;
    VectorXd x0 = d->x0;
    int num_waypoints = d->num_waypoints;
    double velocity_bound = d->velocity_bound;
    if (i == 0) {
        return std::abs((x[3] - x0(1))/(x[1]/num_waypoints)) - velocity_bound;
    } else {
        return std::abs((x[2*(i+1)+1] - x[1+2*i])/(x[1]/num_waypoints)) - velocity_bound;
    }
}

int main(int argc, char const *argv[]) {
    int num_waypoints = 20;
    double steering_velocity_bound = 20;
    double bogie_velocity_bound = 20;
    VectorXd x0 = g1(t1);

    opt optimization = opt(LN_COBYLA, 2+num_waypoints*2);
    optimization.set_xtol_rel(1e-2);
    optimization.set_maxtime(1);

    // Add velocity constraints
    std::vector<VelocityConstraintData> velocityConstraints(num_waypoints);
    for (int i = 0; i < num_waypoints; i++) {
        VelocityConstraintData data = {i, x0, num_waypoints, steering_velocity_bound};
        velocityConstraints[i] = data;
        optimization.add_inequality_constraint(steering_velocity_constraint, &velocityConstraints[i], 1e-8);
        optimization.add_inequality_constraint(bogie_velocity_constraint, &velocityConstraints[i], 1e-8);
    }

    // make initial guess
    double tt_guess = 5;
    std::vector<double> guess = {M_PI, tt_guess};
    VectorXd xf = g2(M_PI);
    VectorXd linear_step = (xf - x0)/num_waypoints;
    for (int i = 1; i < num_waypoints+1; i++) {
        VectorXd waypoint = x0 + linear_step*i;
        guess.push_back(waypoint(0));
        guess.push_back(waypoint(1));
    }

    optimization.set_min_objective(cost_func, NULL);

    double obj_value;
    try {
        result res = optimization.optimize(guess, obj_value);
        std::printf("Optimized Vector {%f, %f}, with objective value %f \n", guess[0], guess[1], obj_value);
    } catch (std::exception &e) {
        std::printf("Optimization failed: %s \n", e.what());
        std::printf("Optimized Vector {%f, %f}, with objective value %f \n", guess[0], guess[1], obj_value);
    }

    std::vector<VectorXd> waypoints(num_waypoints+1, VectorXd(2));
    waypoints[0] = x0;
    for (int i = 1; i < num_waypoints+1; i++) {
        waypoints[i] << guess[2*i], guess[1+2*i];
    }

    Bezier::Spline<double> spline = Bezier::Spline<double>(waypoints, guess[1], t1);
    //plot output
    std::vector<double> xVals, g1Vals, g2Vals, splineVals;
    int num_points = 100;
    double time_delta = guess[1]/(num_points);
    for (int i = 0; i < num_points; i++) {
        double x = t1 + time_delta*i;
        xVals.push_back(x);
        g1Vals.push_back(g1(x)(0));
        g2Vals.push_back(g2(x - t1 - guess[1] + guess[0])(0));
        splineVals.push_back(spline.evaluate(x)(0));
    }

    plt::plot(xVals, g1Vals, "r");
    plt::plot(xVals, g2Vals, "b");
    plt::plot(xVals, splineVals, "g");
    plt::show();

    return 0;
}
