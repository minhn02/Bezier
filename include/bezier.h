#pragma once

#include <iostream>
#include <Eigen/Dense>
#include <vector>
#include <functional>
#include <../src/util.h>
#include <cmath>
#include <assert.h>
#include <map>

using namespace Eigen;

namespace Bezier {

    /**
     * @brief Container for a Bezier Curve
     * evaluation -- https://en.wikipedia.org/wiki/B%C3%A9zier_curve
     * derivative -- https://pages.mtu.edu/~shene/COURSES/cs3621/NOTES/spline/Bezier/bezier-der.html
    */
    class Curve {
        public:

        /**
         * @brief Instantiate a n'th order Bezier Curve defined by n control points over [0, 1]
         * @param pointList List of points that define a bezier curve
        */
        Curve(std::vector<VectorXd> pointList);

        /**
         * @brief Instantiate a n'th order Bezier Curve defined by n control points over [0,T]
         * @param pointList List of points that define a bezier curve
         * @param T the interval [0, T] over which the curve is defined
        */
        Curve(std::vector<VectorXd> pointList, double T);

        /**
         * Evaluates the curve at t
        */
        VectorXd evaluate(double t);

        /**
         * Calculates first derivative of the curve at t
        */
        VectorXd dEvaluate(double t);
        
        /**
         * Calculates the n'th derivative of the curve at t
        */
        VectorXd dEvaluate(int n, double t);

        /**
         * @brief 
         * @param pointNumber point number to change, where 0 is the (1-t) term and n is the t term
         * @param point replace point at pointNumber to this value 
        */
        void setCurvePoint(int pointNumber, VectorXd point);
        
        private:
        /**
         * Returns a new curve that is the derivative of the current curve
        */
        Curve generateDerivativeCurve();

        private:
        std::vector<VectorXd> pointList_;
        int dim_;
        int order_;
        double T_;
        std::vector<double> coefficients_;
        std::map<int, Curve> dMap_;
    };

    class Spline {

        /**
         * Interpolates between points which are over [0, T] with an order degree bezier curve
        */
        Spline(std::vector<VectorXd> points, double T, int order);

        /**
         * Splits func which is defined over [0, T] into N samples, and interpolates between them with an order degree bezier curve
        */
        Spline(std::function<double(double)> func, double T, int N, int order);

        VectorXd evaluate(double t);
        VectorXd dEvaluate(double t);
        VectorXd dEvaluate(int n, double t);
    };
}
