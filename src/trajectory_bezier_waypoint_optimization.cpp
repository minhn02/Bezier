#include "bezier.h"
#include "nlopt.hpp"
#include "nlopt.h"
#include <vector>
#include <cmath>
#include <iostream>
#include "matplotlibcpp.h"
#include "rover_trajectory.h"

using namespace nlopt;
using namespace Eigen;
namespace plt = matplotlibcpp;

// define gaits
double period = 2*M_PI;
auto g1 = [](double t) {VectorXd val(2); val << std::sin(t), std::cos(t); return val;};
auto d_g1 = [](double t) {VectorXd val(2); val << std::cos(t), -std::sin(t); return val;};

auto g2 = [](double t) {VectorXd val(2); val << std::cos(4*t), std::sin(4*t); return val;};
auto d_g2 = [](double t) {VectorXd val(2); val << 4*-std::sin(4*t), 4*std::cos(4*t); return val;};

Bezier::Spline<double> g1_traj = RoverTrajectory::translate_to_cartesian_gait(g1, 1000, period, 0, {0, 0, 0});
Bezier::Spline<double> g2_traj = RoverTrajectory::translate_to_cartesian_gait(g2, 1000, period, 0, {0, 0, 0});

Vector3d g2_traj_displacement;


double t1 = M_PI;
int n_timesteps = 20;
double velocity_constraint = 5;

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
    obj_value += 500*final_diff;

    // penalize points from bezier trajectory
    VectorXd P0 = x0;
    VectorXd P3 = xf;
    VectorXd P1 = P0 + d_g1(t1)*x[1]/3;
    VectorXd P2 = P3 - d_g2(x[0])*x[1]/3;
    Bezier::Curve<double> curve({P0, P1, P2, P3}, x[1], t1);
    double time_delta = x[1]/n_timesteps;

    for (int i = 1; i < n_timesteps+1; i++) {
        Eigen::VectorXd waypoint(2); waypoint << x[2*i], x[2*i+1];
        Eigen::VectorXd bezier_waypoint = curve.evaluate(t1 + time_delta*i);
        obj_value += (bezier_waypoint - waypoint).norm();
    }

    //calculate displacement caused by taking trajectory
    std::vector<VectorXd> waypoints(20+1, VectorXd(2));
    waypoints[0] = x0;
    for (int i = 1; i < 20+1; i++) {
        waypoints[i] << x[2*i], x[1+2*i];
    }
    Bezier::Spline<double> spline = Bezier::Spline<double>(waypoints, x[1], t1);

    size_t num_trajectory_points = 100;
    std::vector<double> beta_list(num_trajectory_points);
    double step = x[1]/num_trajectory_points;
    for (size_t i = 0; i < num_trajectory_points; i++) {
        beta_list[i] = spline.evaluate(t1 + step*i)(0);
    }

    VectorXd initial_vec = g1_traj.evaluate(t1);
    std::vector<double> displacement = RoverTrajectory::calculate_rover_displacement(beta_list, x[1]/(double)num_trajectory_points, {initial_vec(0), initial_vec(1), initial_vec(2)}, x0(0));
    Vector3d displacement_vec; displacement_vec << displacement[0], displacement[1], displacement[2];
    VectorXd actual_displacement = g2_traj.evaluate(x[0]);
    std::cout << "rover displacement: [" << displacement[0] << ", " << displacement[1] << ", " << displacement[2] << "]" << std::endl;
    std::cout << "gait displacement: [" << actual_displacement[0] << ", " << actual_displacement[1] << ", " << actual_displacement[2] << "]" << std::endl;

    // find closest point to rover displacement
    Vector3d closest_point = RoverTrajectory::find_closest_waypoint(displacement_vec, actual_displacement);
    std::cout << "closest_point: [" << closest_point[0] << ", " << closest_point[1] << ", " << closest_point[2] << "]" << std::endl;
    g2_traj_displacement = closest_point - actual_displacement;
    double diff = (closest_point - displacement_vec).norm();
    obj_value += diff;
    std::cout << "displacement_norm: " << (closest_point - displacement_vec).norm() << std::endl;

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
    std::vector<VelocityConstraintData> velocityConstraints(num_waypoints*2);
    for (int i = 0; i < num_waypoints; i+=2) {
        VelocityConstraintData data = {i, x0, num_waypoints, steering_velocity_bound};
        velocityConstraints[i] = data;
        optimization.add_inequality_constraint(steering_velocity_constraint, &velocityConstraints[i], 1e-8);
        data = {i, x0, num_waypoints, bogie_velocity_bound};
        velocityConstraints[i+1] = data;
        optimization.add_inequality_constraint(bogie_velocity_constraint, &velocityConstraints[i], 1e-8);
    }

    // make initial guess
    double tt_guess = 5;
    std::vector<double> guess = {M_PI, tt_guess};
    VectorXd xf = g2(M_PI);

    VectorXd P0 = x0;
    VectorXd P3 = xf;
    VectorXd P1 = P0 + d_g1(t1)*tt_guess/3;
    VectorXd P2 = P3 - d_g2(M_PI)*tt_guess/3;
    Bezier::Curve<double> curve({P0, P1, P2, P3}, tt_guess, t1);

    double guess_time_delta = tt_guess/num_waypoints;
    for (int i = 1; i < num_waypoints+1; i++) {
        VectorXd waypoint = curve.evaluate(t1 + guess_time_delta*i);
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
    int num_points = 1000;
    double time_delta = guess[1]/(num_points - 1);
    for (int i = 0; i < num_points; i++) {
        double x = t1 + time_delta*i;
        xVals.push_back(x);
        g1Vals.push_back(g1(x)(0));
        g2Vals.push_back(g2(x - t1 - guess[1] + guess[0])(0));
        splineVals.push_back(spline.evaluate(x)(0));
    }

    plt::figure();
    plt::named_plot("gait1 joint", xVals, g1Vals, "r");
    plt::named_plot("gait2 joint", xVals, g2Vals, "b");
    plt::named_plot("generated transition", xVals, splineVals, "g");
    plt::legend();
    plt::title("Joint Space Transition");
    plt::show();

    // plot trajectory

    // generate transition trajectory
    VectorXd initial_vector = g1_traj.evaluate(t1);

    std::vector<double> transition_trajectory_vals_x;
    std::vector<double> transition_trajectory_vals_y;

    std::vector<double> g1_trajectory_vals_x;
    std::vector<double> g1_trajectory_vals_y;

    std::vector<double> g2_trajectory_vals_x;
    std::vector<double> g2_trajectory_vals_y;

    std::vector<std::vector<double>> transition_trajectory = RoverTrajectory::transform_joint_movement(splineVals, time_delta, {initial_vector(0), initial_vector(1), initial_vector(2)}, g1(t1)(0));

    double trajectory_time_step = period/(num_points - 1);
    for (int i = 0; i < num_points; i++) {
        transition_trajectory_vals_x.push_back(transition_trajectory[i][0]);
        transition_trajectory_vals_y.push_back(transition_trajectory[i][1]);

        // if (trajectory_time_step*i < t1) {
            g1_trajectory_vals_x.push_back(g1_traj.evaluate(trajectory_time_step*i)(0));
            g1_trajectory_vals_y.push_back(g1_traj.evaluate(trajectory_time_step*i)(1));
        // }

        if (trajectory_time_step*i > guess[0]) {
            g2_trajectory_vals_x.push_back(g2_traj.evaluate(trajectory_time_step*i)(0) + g2_traj_displacement(0));
            g2_trajectory_vals_y.push_back(g2_traj.evaluate(trajectory_time_step*i)(1) + g2_traj_displacement(1));
        }
    }

    plt::figure();
    plt::named_plot("generated transition trajectory", transition_trajectory_vals_x, transition_trajectory_vals_y);
    plt::named_plot("gait1 trajectory", g1_trajectory_vals_x, g1_trajectory_vals_y);
    plt::named_plot("gait2 trajectory", g2_trajectory_vals_x, g2_trajectory_vals_y);
    plt::title("Cartesian Trajectory Transition");
    plt::legend();
    plt::show();

    return 0;
}