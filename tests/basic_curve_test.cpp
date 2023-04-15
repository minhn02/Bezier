#include <cassert>
#include "bezier.h"

/**
 * Points evaluated using python script bezier_calculator.py included in test directory
*/

using namespace Eigen;

void test_cubic_curve_evaluation() {    
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

void test_quintic_curve_evaluation() {
    // initialize vectors
    VectorXd P0(2); P0 << 15, 43;
    VectorXd P1(2); P1 << 1, -5;
    VectorXd P2(2); P2 << 10, -10;
    VectorXd P3(2); P3 << 50, 32;
    VectorXd P4(2); P4 << -79, 23;
    VectorXd P5(2); P5 << 5, -8;

    Bezier::Curve curve({P0, P1, P2, P3, P4, P5});

    // check output vectors, round to 3 decimal places; points evaluated using python script, see top of file
    VectorXd curve0(2); curve0 << 15, 43;
    VectorXd curve01(2); curve01 << 10.284, 23.291;
    VectorXd curve02(2); curve02 << 9.429, 11.777;
    VectorXd curve03(2); curve03 << 10.356, 7.205;
    VectorXd curve04(2); curve04 << 10.386, 7.649;
    VectorXd curve05(2); curve05 << 7.188, 10.781;
    VectorXd curve06(2); curve06 << -0.274, 14.151;
    VectorXd curve07(2); curve07 << -10.789, 15.457;
    VectorXd curve08(2); curve08 << -19.957, 12.823;
    VectorXd curve09(2); curve09 << -19.237, 5.071;
    VectorXd curve1(2); curve1 << 5, -8;

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
    test_cubic_curve_evaluation();
    test_quintic_curve_evaluation();
}
