#pragma once

#include "bezier.h"
#include <algorithm>
#include <iostream>
#include <cmath>
#include <vector>

#include <boost/numeric/odeint.hpp>
#include <Eigen/Dense>


namespace RoverTrajectory {
    /***
        * @brief Returns the displacement of a rover over time
        * @param func a function that returns the velocity of the rover at time t
        * @param num_points the number of points to sample the function at
        * @param t0 the start time of the trajectory
        * @param duration the duration of the trajectory
        * @param initial_position the initial position of the rover in [X, Y, Theta] space
        * @return a vector of [X, Y, Theta] at the end of the trajectory
    */
    std::vector<double> calculate_rover_displacement(std::vector<double> beta_list, double dt, std::vector<double> initial_position, double initial_beta);

    /***
        * @brief Returns the displacement of a rover following a gait trajectory
        * @param func a gait which returns the joint positions at time t
        * @param num_points the number of points to sample the function at
        * @param t0 the start time of the trajectory
        * @param duration the duration of the trajectory
        * @param initial_position the initial position of the rover in [X, Y, Theta] space
        * @return a vector of [X, Y, Theta] at the end of the trajectory
    */
    Bezier::Spline<double> translate_to_cartesian_gait(VectorXd (*func)(double), size_t num_points, double period, double t0, std::vector<double> initial_position);
};
