#include "bezier.h"

namespace Bezier {

    Curve::Curve() {}

    Curve::Curve(std::vector<VectorXd> pointList) : Curve(pointList, 1, 0) {}

    Curve::Curve(std::vector<VectorXd> pointList, double T): Curve(pointList, T, 0) {}

    Curve::Curve(std::vector<VectorXd> pointList, double T, double t0) {
        assert(pointList.size() > 0);

        pointList_ = pointList;
        order_ = pointList.size();
        dim_ = pointList[0].size();
        T_ = T;
        t0_ = t0;
        coefficients_ = Util::generateBinomialCoefficients(order_ - 1);
    }

    VectorXd Curve::evaluate(double t) {
        VectorXd runningSum = VectorXd::Zero(dim_);
        for (int i = 0; i < order_; i++) {
            runningSum += coefficients_[i] *
                          std::pow(((T_ - (t-t0_))/T_), (order_ - 1 - i)) *
                          std::pow((t - t0_)/(T_), i) *
                          pointList_[i];
        }
        return runningSum;
    }

    VectorXd Curve::dEvaluate(double t) {
        return dEvaluate(1, t);
    }

    //TODO cache lower derivatives too (have it return a curve instead of just a value?)
    VectorXd Curve::dEvaluate(int n, double t) {
        assert(order_ >= n);

        if (n == 0) {
            return evaluate(t);
        }

        Curve derivative = this->generateDerivativeCurve();
        if (n == 1) {
            return derivative.evaluate(t);
        } else {
            return derivative.dEvaluate(n - 1, t);
        }
    }

    void Curve::setCurvePoint(int pointNumber, VectorXd point) {
        assert(pointNumber <= order_);

        pointList_[pointNumber] = point;
    }

    Curve Curve::generateDerivativeCurve() {
        if (!dMap_.count(1)) {
            std::vector<VectorXd> pointList(order_ - 1);
            for (int i = 0; i < order_ - 1; i++) {
                VectorXd point = (order_ - 1) * (pointList_[i+1] - pointList_[i]);
                pointList[i] = point;
            }
            dMap_.insert({1, Curve(pointList, T_, t0_)});
        }

        return dMap_.at(1);
    }

    Spline::Spline(std::vector<Curve> curves, double T, double t0) {
        curves_ = curves;
        T_ = T;
        t0_ = t0;
    }

    Spline::Spline(std::vector<VectorXd> points, double T, double t0) {
        T_ = T;
        t0_ = t0;
        initializePointList(points, T, t0);
    }
    
    Spline::Spline(std::function<VectorXd(double)> func, int N, double T, double t0) {
        T_ = T;
        t0_ = t0;
        std::vector<VectorXd> points(N);
        double delta = T/((double)N);
        int index = 0;
        for (double i = t0_; i <= t0_ + T; i += delta) {
            points[index] = func(i);
            index += 1;
        }
        assert(index == N);
        initializePointList(points, T, t0);
    }

    void Spline::initializePointList(std::vector<VectorXd> points, double T, double t0) {
        T_ = T;
        t0_ = t0;
        int n = points.size() - 1;
        std::vector<double> uk(n+1);
        std::vector<VectorXd> tk(n+1);
        curves_.resize(n);
        int dim = points[0].size();

        // uk (time delta between points) from (8.13) pg 366 Trajectory Planning For Automatic Machines and Robots
        // tk (tangent vectors at each point) from (8.49) pg 400 of Trajectory Planning For Automatic Machines and Robots
        // calculate u_k according to equally spaced distribution
        double d = 0;
        for (int k = 1; k < points.size(); k++) {
            d += (points[k] - points[k-1]).norm();
        }
        uk[0] = 0;
        for (int i = 1; i < n + 1; i++) {
            uk[i] = ((double)i/(double)n)*T_;
        }

        // calculate tk
        for (int k = 1; k < n; k++) {
            Eigen::VectorXd delta_q_k = points[k] - points[k-1];
            Eigen::VectorXd delta_q_k_1 = points[k+1] - points[k];
            double delta_u_k = uk[k] - uk[k-1];
            double delta_u_k_1 = uk[k+1] - uk[k];
            double alpha_k = delta_u_k / (delta_u_k + delta_u_k_1);

            tk[k] = (1 - alpha_k)*(delta_q_k/delta_u_k) + alpha_k*(delta_q_k_1/delta_u_k_1);
        }
        // calculate endpoints
        tk[0] = 2 * ((points[1] - points[0]) / (uk[1] - uk[0])) - tk[1];
        tk[n] = 2 * ((points[n] - points[n-1]) / (uk[n] - uk[n-1])) - tk[n-1];

        // bezier curves calculated from (8.51) pg 401 Trajectory Planning For Automatic Machines and Robots
        // yields n-1 curves where len(points) = n
        for (int i = 0; i < n; i++) {
            // compute alpha
            double alpha = calculateTangentMagnitude(points[i], points[i+1], tk[i], tk[i+1]);

            // consruct curve
            VectorXd P0 = points[i];
            VectorXd P3 = points[i+1];
            VectorXd P1 = P0 + (1/3)*alpha*tk[i];
            VectorXd P2 = P3 - (1/3)*alpha*tk[i+1];
            curves_[i] = Curve({P0, P1, P2, P3}, uk[i+1] - uk[i], uk[i]);
        }
    }

    /**
     * (8.52) pg 402 of trajectory textbook
    */
    double Spline::calculateTangentMagnitude(VectorXd p0, VectorXd p3, VectorXd t0, VectorXd t3) {
        double a = 16 - std::pow((t0 + t3).norm(), 2);
        double b = 12 * (p3 - p0).transpose().dot(t0 + t3);
        double c = -36 * std::pow((p3 - p0).norm(), 2);

        double discriminant = std::pow(b, 2) - 4 * a * c;

        if (discriminant < 0) {
            return std::numeric_limits<double>::quiet_NaN();
        }

        double alpha1 = (-b + std::sqrt(discriminant)) / (2 * a);
        double alpha2 = (-b - std::sqrt(discriminant)) / (2 * a);

        // Return the positive solution for alpha
        return alpha1 > 0 ? alpha1 : alpha2;     
    }


    VectorXd Spline::evaluate(double t) {
        // calculate curve in spline to evaluate
        //TODO undefined if t is outside [0, T]
        t = t - t0_;
        int n = curves_.size() - 1;
        int index = std::round((t / T_) * (n));

        return curves_[index].evaluate(t);
    }

    VectorXd Spline::dEvaluate(double t) {
        return dEvaluate(1, t);
    }

    VectorXd Spline::dEvaluate(int n, double t) {
        // calculate curve in spline to evaluate
        t = t - t0_;
        int len = curves_.size() - 1;
        int index = std::round((t / T_) * (len));

        return curves_[index].dEvaluate(n, t);
    }
}