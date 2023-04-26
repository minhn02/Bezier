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
        initializePointList(points, T, t0);
    }
    
    Spline::Spline(std::function<VectorXd(double)> func, int N, double T, double t0) {
        std::vector<VectorXd> points(N);
        double delta = T/(N-1);
        int index = 0;
        for (double i = 0; i <= T; i += delta) {
            points[index] = func(i);
            index += 1;
        }
        initializePointList(points, T, t0);
    }

    void Spline::initializePointList(std::vector<VectorXd> points, double T, double t0) {
        T_ = T;
        t0_ = t0;
        int n = points.size() - 1;
        std::vector<double> uk(n+1);
        std::vector<VectorXd> tk(n+1);
        curves_.resize(n-1);
        int dim = points[0].size();

        // uk (time delta between points) from (8.13) pg 366 Trajectory Planning For Automatic Machines and Robots
        // tk (tangent vectors at each point) from (8.49) pg 400 of Trajectory Planning For Automatic Machines and Robots
        // calculate u_k according to cord length distribution
        double d = 0;
        for (int k = 1; k < n; k++) {
            d += (points[k] - points[k-1]).norm();
        }
        uk[0] = 0;
        uk[n] = T_;
        for (int i = 1; i < n; i++) {
            uk[i] = uk[i-1] + ((points[i] - points[i-1]).norm() / d) * T_;
        }

        std::vector<VectorXd> deltak(n+1);
        std::vector<double> alphak(n);
        // calculate deltak and alphak
        for (int k = 1; k < n+1; k++) {
            deltak[k] = (points[k] - points[k-1]).array() / (uk[k] - uk[k-1]);
            // alphak only defined for 1 <= k <= n-1
            if (k != n) {
                alphak[k] = (uk[k] - uk[k-1]) / ((uk[k] - uk[k-1]) + (uk[k+1] - uk[k]));
            }
        }

        // calculate tk
        for (int k = 1; k < n; k++) {
            tk[k] = (1 - alphak[k])*deltak[k] + (alphak[k]*deltak[k+1]);
        }
        // calculate endpoints
        tk[0] = 2*deltak[1] - tk[1];
        tk[n] = 2*deltak[n] - tk[n-1];

        // bezier curves calculated from (8.51) pg 401 Trajectory Planning For Automatic Machines and Robots
        // yields n-1 curves where len(points) = n
        for (int i = 0; i < n-1; i++) {
            // compute alpha
            VectorXd a = 16 - (tk[i] + tk[i+1]).cwiseAbs2().array();
            VectorXd b = 12 * (points[i+1] - points[i]) * (tk[i] + tk[i+1]);
            VectorXd c = -36 * (points[i+1] - points[i]).cwiseAbs2();
            VectorXd delta = b * b - 4 * a * c;
            VectorXd alpha = (-1*b.array() - delta.array().sqrt()) / (2 * a.array());

            // consruct curve
            VectorXd P0 = points[i];
            VectorXd P3 = points[i+1];
            VectorXd P1 = P0 + (1/3)*alpha*tk[i];
            VectorXd P2 = P3 - (1/3)*alpha*tk[i+1];
            curves_[i] = Curve({P0, P1, P2, P3}, uk[i]);
        }
    }

    VectorXd Spline::evaluate(double t) {
        // calculate curve in spline to evaluate
        int n = curves_.size();
        int index = std::floor((t / T_) * (n-1));

        return curves_[index].evaluate(t);
    }

    VectorXd Spline::dEvaluate(double t) {
        return dEvaluate(1, t);
    }

    VectorXd Spline::dEvaluate(int n, double t) {
        // calculate curve in spline to evaluate
        int index = std::floor(T_ / t);

        return curves_[index].dEvaluate(n, t);
    }
}