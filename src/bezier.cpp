#include "bezier.h"

namespace Bezier {

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

    Spline::Spline(std::vector<Curve> curves, double T) {
        curves_ = curves;
        T_ = T;
    }

    Spline::Spline(std::vector<VectorXd> points, double T) {
        T_ = T;
        int n = points.size();

        std::vector<Curve> curveList;

        // u_k (time delta between points) from (8.13) pg 514 Trajectory Planning For Automatic Machines and Robots
        // (8.49) pg 400 of Trajectory Planning For Automatic Machines and Robots
        double d = 0;
        for (int k = 1; k < n; k++) {
            d += (points[k] - points[k-1]).norm();
        }

        double u_k_1 = 0; 
        for (int i = 1; i < n; i++) {
            double u_k = u_k_1 + (points[i] - points[i-1]).norm() / d;
            u_k_1 = u_k;

            VectorXd P0 = points[i-1];
            VectorXd P3 = points[i];

        }
    }

    VectorXd Spline::evaluate(double t) {
        // calculate curve in spline to evaluate
        int index = std::floor(T_ / t);

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