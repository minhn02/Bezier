#pragma once

#include <iostream>
#include <Eigen/Dense>
#include <vector>
#include <functional>
#include <cmath>
#include <assert.h>
#include <map>

using namespace Eigen;

namespace Bezier {

    /**
    * @brief Returns binomial coefficients [(n 0), (n 1), ..., (n n)]
    */
    inline std::vector<double> generateBinomialCoefficients(int n) {
        std::vector<double> row = {1};
        for (int k=0; k < n; k++) {
            row.push_back(row[k] * (n-k) / (k+1));
        }
        return row;
    }

    /**
     * @brief Container for a Bezier Curve
     * evaluation -- https://en.wikipedia.org/wiki/B%C3%A9zier_curve
     * derivative -- https://pages.mtu.edu/~shene/COURSES/cs3621/NOTES/spline/Bezier/bezier-der.html
    */
   template <typename X>
    class Curve {
        public:

        /**
         * @brief creates an empty Bezier Curve container
        */
        Curve() {}

        /**
         * @brief Instantiate a n'th order Bezier Curve defined by n control points over [t0,T+t0]
         * @param pointList List of points that define a bezier curve
         * @param T the interval [t0, T+t0] over which the curve is defined
         * @param t0 the start time of the curve
        */
        Curve(std::vector<VectorXd> pointList, X T=1, X t0=0) {
            assert(pointList.size() > 0);

            pointList_ = pointList;
            order_ = pointList.size();
            dim_ = pointList[0].size();
            T_ = T;
            t0_ = t0;
            coefficients_ = generateBinomialCoefficients(order_ - 1);
        }

        /**
         * Evaluates the curve at t
        */
        VectorXd evaluate(X t) {
            VectorXd runningSum = VectorXd::Zero(dim_);
            for (int i = 0; i < order_; i++) {
                runningSum += coefficients_[i] *
                              std::pow(1 - ((t-t0_)/(double)T_), (order_ - 1 - i)) *
                              std::pow((t - t0_)/(double)T_, i) *
                              pointList_[i];
            }
            return runningSum;
        }

        /**
         * Calculates first derivative of the curve at t
        */
        VectorXd dEvaluate(X t) {
            return dEvaluate(1, t);
        }
        
        /**
         * Calculates the n'th derivative of the curve at t
        */
        VectorXd dEvaluate(int n, X t) {
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

        /**
         * @brief 
         * @param pointNumber point number to change, where 0 is the (1-t) term and n is the t term
         * @param point replace point at pointNumber to this value 
        */
        void setCurvePoint(int pointNumber, VectorXd point) {
            assert(pointNumber >= 0 && pointNumber < order_);
            assert(point.size() == dim_);
            pointList_[pointNumber] = point;
        }

        X getT() {
            return T_;
        }
        X getT0() {
            return t0_;
        }
        std::vector<VectorXd> getPoints() {
            return pointList_;
        }

        bool isFinished(X t) {
            return t > t0_ + T_;
        }
        
        private:
        /**
         * Returns a new curve that is the derivative of the current curve
        */
        Curve generateDerivativeCurve() {
            std::vector<VectorXd> derivativePointList;
            for (int i = 0; i < order_ - 1; i++) {
                derivativePointList.push_back((1.0/T_) *(order_ - 1) * (pointList_[i+1] - pointList_[i]));
            }
            return Curve(derivativePointList, T_, t0_);
        }

        std::vector<VectorXd> pointList_;
        int dim_;
        int order_;
        X T_;
        X t0_;
        std::vector<double> coefficients_;
        std::map<int, Curve> dMap_;
    };

    template <typename X>
    class Spline {

        public:

        Spline() {
        }

        /**
         * Curves that form a continous spline over [t0, T+t0]
         * Each curve is spaced so it takes up T/len(curves) time
        */
        Spline(std::vector<Curve<X>> curves, X T=1, X t0=0) {
            assert(curves.size() > 0);
            curves_ = curves;
            T_ = T;
            t0_ = t0;
        }

        /**
         * Interpolates between points which are over [0, T] with a 3rd degree bezier curve
        */
        Spline(std::vector<VectorXd> points, X T=1, X t0=0) {
            T_ = T;
            t0_ = t0;
            initializePointList(points, T, t0);
        }

        /**
         * Splits func which is defined over [0, T] into N samples, and interpolates between them with a 3rd degree bezier curve defined over [t0, T+t0]
        */
        Spline(std::function<VectorXd(X)> func, int N, X T=1, X t0=0) {
            T_ = T;
            t0_ = t0;
            std::vector<VectorXd> points(N);
            double delta = T/((double)(N-1));
            double t = t0_;
            for (int i = 0; i < N; i++) {
                points[i] = func(t);
                t += delta;
            }
            initializePointList(points, T, t0);
        }

        VectorXd evaluate(X t) {
            // calculate curve in spline to evaluate
            //TODO undefined if t is outside [0, T]
            t = t - t0_;
            int n = curves_.size();
            int index = std::floor(((double)t / (double)T_) * (n));

            if (index == n) {
                return curves_[n-1].evaluate(t0_ + T_);
            }

            return curves_[index].evaluate(t);
        }
        VectorXd dEvaluate(X t) {
            return dEvaluate(1, t);
        }
        VectorXd dEvaluate(int n, X t) {
            // calculate curve in spline to evaluate
            t = t - t0_;
            int len = curves_.size();
            int index = std::floor(((double)t / (double)T_) * (len));

            return curves_[index].dEvaluate(n, t);
        }

         X getT() {
            return T_;
        }
        X getT0() {
            return t0_;
        }

        bool isFinished(X t) {
            return t > t0_ + T_;
        }

        private:
        void initializePointList(std::vector<VectorXd> points, X T, X t0) {
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
            // calculate tk at endpoints
            tk[0] = 2 * ((points[1] - points[0]) / (uk[1] - uk[0])) - tk[1];
            tk[n] = 2 * ((points[n] - points[n-1]) / (uk[n] - uk[n-1])) - tk[n-1];

            // bezier curves calculated from (8.51) pg 401 Trajectory Planning For Automatic Machines and Robots
            // yields n-1 curves where len(points) = n
            for (int i = 0; i < n; i++) {
                double delta = uk[i+1] - uk[i];
                // consruct curve
                VectorXd P0 = points[i];
                VectorXd P3 = points[i+1];
                VectorXd P1 = P0 + ((double)1/(double)3)*delta*tk[i];
                VectorXd P2 = P3 - ((double)1/(double)3)*delta*tk[i+1];
                curves_[i] = Curve<X>({P0, P1, P2, P3}, delta, uk[i]);
            }
        }

        std::vector<Curve<X>> curves_;
        X T_;
        X t0_;
    };
}
