#include <cassert>
#include "bezier.h"

/**
 * Points evaluated using python script bezier_calculator.py included in test directory
*/

using namespace Eigen;

void test_d_cubic_curve_evaluation() {    
    // initialize vectors
    VectorXd P0(2); P0 << 10, 5;
    VectorXd P1(2); P1 << 8, 2;
    VectorXd P2(2); P2 << 6, 3;
    VectorXd P3(2); P3 << 1, 20;

    Bezier::Curve curve({P0, P1, P2, P3});

    // check output vectors, round to 3 decimal places; points evaluated using python script, see top of file
    VectorXd curve0(2); curve0 << -6, -9;
    VectorXd curve01(2); curve01 << -6.09, -6.24;
    VectorXd curve02(2); curve02 << -6.36, -2.76;
    VectorXd curve03(2); curve03 << -6.81, 1.44;
    VectorXd curve04(2); curve04 << -7.44, 6.36;
    VectorXd curve05(2); curve05 << -8.25, 12;
    VectorXd curve06(2); curve06 << -9.24, 18.36;
    VectorXd curve07(2); curve07 << -10.41, 25.44;
    VectorXd curve08(2); curve08 << -11.76, 33.24;
    VectorXd curve09(2); curve09 << -13.29, 41.76;
    VectorXd curve1(2); curve1 << -15, 51;

    assert((curve.dEvaluate(0) - curve0).norm() < 0.01);
    assert((curve.dEvaluate(0.1) - curve01).norm() < 0.01);
    assert((curve.dEvaluate(0.2) - curve02).norm() < 0.01);
    assert((curve.dEvaluate(0.3) - curve03).norm() < 0.01);
    assert((curve.dEvaluate(0.4) - curve04).norm() < 0.01);
    assert((curve.dEvaluate(0.5) - curve05).norm() < 0.01);
    assert((curve.dEvaluate(0.6) - curve06).norm() < 0.01);
    assert((curve.dEvaluate(0.7) - curve07).norm() < 0.01);
    assert((curve.dEvaluate(0.8) - curve08).norm() < 0.01);
    assert((curve.dEvaluate(0.9) - curve09).norm() < 0.01);
    assert((curve.dEvaluate(1) - curve1).norm() < 0.01);
}


//TODO fill out correct values
void test_d_d_cubic_curve_evaluation() {    
    // initialize vectors
    VectorXd P0(2); P0 << 10, 5;
    VectorXd P1(2); P1 << 8, 2;
    VectorXd P2(2); P2 << 6, 3;
    VectorXd P3(2); P3 << 1, 20;

    Bezier::Curve curve({P0, P1, P2, P3});

    // check output vectors, round to 3 decimal places; points evaluated using python script, see top of file
    VectorXd curve0(2); curve0 << 10, 5;
    VectorXd curve01(2); curve01 << 9.397, 4.232;
    VectorXd curve02(2); curve02 << 8.776, 3.776;
    VectorXd curve03(2); curve03 << 8.119, 3.704;
    VectorXd curve04(2); curve04 << 7.408, 4.088;
    VectorXd curve05(2); curve05 << 6.625, 5;
    VectorXd curve06(2); curve06 << 5.752, 6.512;
    VectorXd curve07(2); curve07 << 4.771, 8.696;
    VectorXd curve08(2); curve08 << 3.664, 11.624;
    VectorXd curve09(2); curve09 << 2.413, 15.368;
    VectorXd curve1(2); curve1 << 1, 20;

    assert((curve.evaluate(0) - curve0).norm() < 0.01);
    assert((curve.evaluate(0.1) - curve01).norm() < 0.01);
    assert((curve.evaluate(0.2) - curve02).norm() < 0.01);
    assert((curve.evaluate(0.3) - curve03).norm() < 0.01);
    assert((curve.evaluate(0.4) - curve04).norm() < 0.01);
    assert((curve.evaluate(0.5) - curve05).norm() < 0.01);
    assert((curve.evaluate(0.6) - curve06).norm() < 0.01);
    assert((curve.evaluate(0.7) - curve07).norm() < 0.01);
    assert((curve.evaluate(0.8) - curve08).norm() < 0.01);
    assert((curve.evaluate(0.9) - curve09).norm() < 0.01);
    assert((curve.evaluate(1) - curve1).norm() < 0.01);
}

int main(int argc, char const *argv[]) {
    test_d_cubic_curve_evaluation();
    // test_d_d_cubic_curve_evaluation();
}
