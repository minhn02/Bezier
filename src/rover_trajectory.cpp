#include "rover_trajectory.h"

using namespace Eigen;

// state is [X, Y, Theta]
// X = X0 + cos(Theta)*X'
// Y = Y0 + sin(Theta)*Y'
// Theta = Theta0 + Vr/lx * tan(beta/2) + d_beta/2
typedef Vector3d state_type;
double DEG_TO_RAD = M_PI/180.0;
double RAD_TO_DEG = 180.0/M_PI;
// constants
double lx = 0.200;
double ly = 0.200;
double r = 0.100;

double dbeta_dt = 15.0;

struct push_back_state_and_time
{
    std::vector< state_type >& m_states;
    std::vector< double >& m_times;

    push_back_state_and_time( std::vector< state_type > &states , std::vector< double > &times )
    : m_states( states ) , m_times( times ) { }

    void operator()( const state_type &x , double t )
    {
        m_states.push_back( x );
        m_times.push_back( t );
    }
};

Matrix3d createYawRotationMatrix(double yaw, bool degrees=false) {
    if (degrees) {
        yaw *= DEG_TO_RAD;
    }

    Matrix3d yaw_mat {
        {std::cos(yaw), -std::sin(yaw), 0},
        {std::sin(yaw), std::cos(yaw), 0},
        {0, 0, 1}
    };

    return yaw_mat;
}

class wheel_displacement {

    double Vw_;
    double Vr_;
    uint8_t wheel_;
    double beta_;
    double dbdt_;

    public:
    wheel_displacement( double Vw, double Vr, uint8_t wheel, double beta, double dbdt ) : Vw_(Vw), Vr_(Vr), wheel_(wheel), beta_(beta), dbdt_(dbdt) { }

    void operator() ( const state_type &x , state_type &dxdt , const double /* t */ ) {
        double theta_ = x(2);

        dxdt << 
        std::cos(theta_)*Vw_,
        std::sin(theta_)*Vw_,
        Vr_/lx * std::tan(-beta_/2) + ((wheel_<2) ? -1 : 1)*dbdt_*DEG_TO_RAD/2;
    }
};

/**
 * @brief based on the movement of the steering angle, calculates the displacement of the rover
 * 
 * @param beta_list list of steering joint angles
 * @param dt time between two points in beta_list
 * @return std::vector<double> [delta_x, delta_y]
 */
std::vector<double> calculate_rover_displacement(std::vector<double> beta_list, double dt) {
    double beta_max = *std::max_element(beta_list.begin(), beta_list.end());
    size_t num_points = beta_list.size();

    std::vector<double> dbdt;
    std::vector<std::vector<double>> w(4, std::vector<double>{});
    std::vector<double> V;
    double beta_prev = beta_list[0];

    Matrix3d Rbeta = createYawRotationMatrix(beta_max/2, true);
    
    std::vector<Vector3d> Xw;
    Vector3d wheel1; wheel1 << lx, ly, beta_max/2*DEG_TO_RAD;
    Vector3d wheel2; wheel2 << lx, -ly, beta_max/2*DEG_TO_RAD;
    Vector3d wheel3; wheel3 << -lx, ly, -beta_max/2*DEG_TO_RAD;
    Vector3d wheel4; wheel4 << -lx, -ly, -beta_max/2*DEG_TO_RAD;

    Xw.push_back(Rbeta * wheel1);
    Xw.push_back(Rbeta * wheel2);
    Xw.push_back(Rbeta.transpose() * wheel3);
    Xw.push_back(Rbeta.transpose() * wheel4);

    for (size_t i = 0; i < beta_list.size(); i++) {
        double beta = beta_list[i]*DEG_TO_RAD;
        short direction;

        if (beta >= beta_prev) {
            direction = 1;
        } else {
            direction = -1;
        }
        dbdt.push_back(direction*dbeta_dt);

        std::vector<double> vs(4);
        std::vector<double> vd(4);
        for (short wheel = 0; wheel < 4; wheel++) {
            vs[wheel] = (((wheel/2) ? -1: 1) * ((-lx)*std::tan(beta/2) + ((wheel%2) ? -1 : 1)*ly)*direction*dbeta_dt*DEG_TO_RAD/2);
            vd[wheel] = (((wheel%2) ? -1 : 1)*ly/lx * std::tan(beta/2));
        }

        double rover_velocity = 0;
        for (short wheel = 0; wheel < 4; wheel++) {
            rover_velocity = std::max(rover_velocity, -vs[wheel]/(1+vd[wheel]));
        }
        V.push_back(rover_velocity);

        for (short wheel = 0; wheel < 4; wheel++) {
            w[wheel].push_back((V.back()*(1+vd[wheel]) + vs[wheel])/r);

            wheel_displacement wheel_disp(r*w[wheel][i], V[i], wheel, beta, dbdt[i]);

            std::vector<state_type> x_vec;
            std::vector<double> times;
            size_t steps = boost::numeric::odeint::integrate(wheel_disp, Xw[wheel], 0.0, dt, dt/100, push_back_state_and_time(x_vec, times));
        }

        beta_prev = beta;
    }

    double delta_x, delta_y = 0;
    for (uint8_t wheel = 0; wheel < 4; wheel++) {
        delta_x += Xw[wheel](0);
        delta_y += Xw[wheel](1);
    }

    return {(delta_x/4)*1e3, (delta_y/4)*1e3};

}

std::vector<double> calculate_rover_displacement(double (*func)(double), size_t num_points, double t0, double duration) {
    // interpolate func with num_points

    std::vector<double> beta_list(num_points);
    double step = duration/(double)num_points;

    for (size_t i = 0; i < num_points; i++) {
        beta_list[i] = func(t0+step*i);
    }

    return calculate_rover_displacement(beta_list, duration/num_points);
}

double sin_gait(double t) {
    return std::sin(t);
}

int main(int argc, char** argv) {
    
    std::vector<double> displacement = calculate_rover_displacement(sin_gait, 1000, 0, 10);
    std::cout << "X: " << displacement[0] << "\tY: " << displacement[1] << std::endl;
    return 0;
}