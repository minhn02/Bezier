#include "rover_trajectory.h"

using namespace Eigen;

// state is [X, Y, Theta]
// X = X0 + cos(Theta)*X'
// Y = Y0 + sin(Theta)*Y'
// Theta = Theta0 + Vr/lx * tan(beta/2) + d_beta/2

namespace RoverTrajectory {
typedef std::vector<double> state_type;
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

        Matrix3d yaw_mat;
        yaw_mat << 
            std::cos(yaw), -std::sin(yaw), 0,
            std::sin(yaw), std::cos(yaw), 0,
            0, 0, 1;

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
        double theta_ = x[2];

        dxdt[0] = std::cos(theta_)*Vw_;
        dxdt[1] = std::sin(theta_)*Vw_;
        dxdt[2] = Vr_/lx * std::tan(-beta_/2) + ((wheel_<2) ? -1 : 1)*dbdt_*DEG_TO_RAD/2;
    }
};

std::vector<std::vector<double>> eval_displacement_points(std::vector<double> beta_list, double dt, std::vector<double> initial_position, double initial_beta) {
    size_t num_points = beta_list.size();
    // flip sign of beta for convention
    initial_beta *= -1;

    std::vector<std::vector<double>> evaluated_displacement_points(num_points);

    std::vector<double> dbdt;
    std::vector<std::vector<double>> w(4, std::vector<double>{});
    std::vector<double> V;
    double beta_prev = beta_list[0];
    double delta_x, delta_y, delta_theta, theta_back, theta_front = 0;

    Matrix3d Rbeta = createYawRotationMatrix(initial_beta/2, true);
    
    std::vector<Vector3d> Xw;
    Vector3d wheel1; wheel1 << lx + initial_position[0], ly + initial_position[1], (initial_beta/2 + initial_position[2])*DEG_TO_RAD;
    Vector3d wheel2; wheel2 << lx + initial_position[0], -ly + initial_position[1], (initial_beta/2 + initial_position[2])*DEG_TO_RAD;
    Vector3d wheel3; wheel3 << -lx + initial_position[0], ly + initial_position[1], -(initial_beta/2 + initial_position[2])*DEG_TO_RAD;
    Vector3d wheel4; wheel4 << -lx + initial_position[0], -ly + initial_position[1], -(initial_beta/2 + initial_position[2])*DEG_TO_RAD;

    Xw.push_back(Rbeta * wheel1);
    Xw.push_back(Rbeta * wheel2);
    Xw.push_back(Rbeta.transpose() * wheel3);
    Xw.push_back(Rbeta.transpose() * wheel4);

    std::vector<std::vector<double>> Xw_vec(4, std::vector<double>{});
    for (size_t i = 0; i < 4; i++) {
        Xw_vec[i].push_back(Xw[i](0));
        Xw_vec[i].push_back(Xw[i](1));
        Xw_vec[i].push_back(Xw[i](2));
    }

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
            size_t steps = boost::numeric::odeint::integrate(wheel_disp, Xw_vec[wheel], 0.0, dt, dt/100, push_back_state_and_time(x_vec, times));
        }

        beta_prev = beta;

        delta_x = 0;
        delta_y = 0;
        delta_theta = 0;
        for (short wheel = 0; wheel < 4; wheel++) {
            delta_x += Xw_vec[wheel][0];
            delta_y += Xw_vec[wheel][1];
        }
        delta_x = (delta_x/4)*1e3;
        delta_y = (delta_y/4)*1e3;

        //find theta perpendicular to the rover's front axle (between wheels 1 and 2)
        double theta_front = std::atan2(Xw_vec[0][1] - Xw_vec[1][1], Xw_vec[1][0] - Xw_vec[0][0]);
        double theta_back = std::atan2(Xw_vec[2][1] - Xw_vec[3][1], Xw_vec[3][0] - Xw_vec[2][0]);
        delta_theta = (theta_front + theta_back)/2;

        evaluated_displacement_points[i] = {delta_x, delta_y, delta_theta};
    }

    return evaluated_displacement_points;
}

std::vector<double> calculate_rover_displacement(std::vector<double> beta_list, double dt, std::vector<double> initial_position, double initial_beta) {
    return eval_displacement_points(beta_list, dt, initial_position, initial_beta).back();
}

//TODO change to passing in a gait and can call .evaluate() to get positions
Bezier::Spline<double> translate_to_cartesian_gait(VectorXd (*func)(double), size_t num_points, double period, double t0, std::vector<double> initial_position) {
    std::vector<double> beta_list(num_points);
    double step = period/(double)num_points;

    for (size_t i = 0; i < num_points; i++) {
        beta_list[i] = func(t0+step*i)(0);
    }

    std::vector<std::vector<double>> displacements = eval_displacement_points(beta_list, step, initial_position, -40.0);
    std::vector<VectorXd> displacement_vectors(displacements.size(), VectorXd(3));

    for (size_t i = 0; i < displacements.size(); i++) {
        displacement_vectors[i] << displacements[i][0], displacements[i][1], displacements[i][2];
    }

    return Bezier::Spline<double>(displacement_vectors, period, t0);
}
};