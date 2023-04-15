#include "bezier.h"

namespace Bezier {

    Curve::Curve(std::vector<VectorXd> pointList) : Curve(pointList, 1) {}

    Curve::Curve(std::vector<VectorXd> pointList, double T) {
        assert(pointList.size() > 0);

        pointList_ = pointList;
        order_ = pointList.size();
        dim_ = pointList[0].size();
        T_ = T;
        coefficients_ = Util::generateBinomialCoefficients(order_ - 1);
    }

    VectorXd Curve::evaluate(double t) {
        VectorXd runningSum = VectorXd::Zero(dim_);
        for (int i = 0; i < order_; i++) {
            runningSum += coefficients_[i] *
                          std::pow(((T_ - t)/T_), (order_ - 1 - i)) *
                          std::pow(t / T_, i) *
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
            //TODO cache recursive calculations of curves to dMap
            return derivative.dEvaluate(n - 1, t);
        }
    }

    void Curve::setCurvePoint(int pointNumber, VectorXd point) {
        assert(pointNumber <= order_);

        pointList_[pointNumber] = point;
    }

    Curve Curve::generateDerivativeCurve() {
        if (!dMap_.count(1)) {
            std::vector<VectorXd> pointList;
            for (int i = 0; i < order_ - 1; i++) {
                VectorXd point = (order_ - 1) * (pointList_[i+1] - pointList_[i]);
                pointList.insert(pointList.begin() + i, point);
            }
            dMap_.insert({1, Curve(pointList)});
        }

        return dMap_.at(1);
    }
}