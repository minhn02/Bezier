#include <algorithm>
#include <iostream>
#include <cmath>
#include <vector>

#include <boost/numeric/odeint.hpp>
#include <Eigen/Dense>

std::vector<double> calculate_rover_displacement(double (*func)(double), size_t num_points, double t0, double duration);