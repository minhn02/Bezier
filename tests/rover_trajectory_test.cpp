#include <cassert>
#include <vector>
#include <cstdio>
#include <chrono>

#include "bezier.h"
#include "rover_trajectory.h"

void test_squirming_trajectory() {
    size_t num_points = 10000;
    double beta_max = 40;
    double dbeta_dt = 15;
    double dt = beta_max*2/num_points/dbeta_dt;

    double step = (beta_max*2)/num_points;
    std::vector<double> beta_list(num_points);
    for (size_t i = 0; i < num_points; i++) {
        beta_list[i] = -beta_max + step*i;
    }

    std::vector<double> displacement = RoverTrajectory::calculate_rover_displacement(beta_list, dt, {0, 0, 0}, -40.0);
    assert(std::abs(displacement[0] - 200.1054) < 0.1);
    assert(std::abs(displacement[1] - 87.726) < 0.1);
    assert(std::abs(displacement[2] - 1.5708) < 0.1);
}

int main(int argc, char const *argv[]) {
    test_squirming_trajectory();
}
