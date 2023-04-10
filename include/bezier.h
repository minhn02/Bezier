#pragma once

#include <iostream>
#include <Eigen/Dense>
#include <vector>

using namespace Eigen;

namespace Bezier {

    /**
     * @brief Container for a Bezier Curve
    */
    class Curve {
        public:

        /**
         * @brief Instantiate a n'th order Bezier Curve defined by n control points
         * @param pointList List of points that define a bezier curve
        */
        Curve(std::vector<VectorXd> pointList);

        VectorXd evaluate(double t);
        VectorXd dEvaluate(double t);
        VectorXd ddEvaluate(double t);

        /**
         * @brief 
         * @param pointNumber point number to change, where 0 is the (1-t) term and n is the t term
         * @param point replace point at pointNumber to this value 
        */
        void setCurvePoint(int pointNumber, VectorXd point);
    };
     
}
