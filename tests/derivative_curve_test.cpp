#include <cassert>
#include "bezier.h"

/**
 * Points evaluated using python script bezier_calculator.py included in test directory
*/

using namespace Eigen;


//TODO redo test cases for derivatives; use chain rule (1/T)*B(t-t0/T)
void test_d_cubic_curve_evaluation() {    
    // initialize vectors
    VectorXd P0(2); P0 << 10, 5;
    VectorXd P1(2); P1 << 8, 2;
    VectorXd P2(2); P2 << 6, 3;
    VectorXd P3(2); P3 << 1, 20;

    Bezier::Curve<double> curve({P0, P1, P2, P3});

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


void test_d_d_cubic_curve_evaluation() {    
    // initialize vectors
    VectorXd P0(2); P0 << 10, 5;
    VectorXd P1(2); P1 << 8, 2;
    VectorXd P2(2); P2 << 6, 3;
    VectorXd P3(2); P3 << 1, 20;

    Bezier::Curve<double> curve({P0, P1, P2, P3});

    // check output vectors, round to 3 decimal places; points evaluated using python script, see top of file
    VectorXd curve0(2); curve0 << 0, 24;
    VectorXd curve01(2); curve01 << -1.8, 31.2;
    VectorXd curve02(2); curve02 << -3.6, 38.4;
    VectorXd curve03(2); curve03 << -5.4, 45.6;
    VectorXd curve04(2); curve04 << -7.2, 52.8;
    VectorXd curve05(2); curve05 << -9, 60;
    VectorXd curve06(2); curve06 << -10.8, 67.2;
    VectorXd curve07(2); curve07 << -12.6, 74.4;
    VectorXd curve08(2); curve08 << -14.4, 81.6;
    VectorXd curve09(2); curve09 << -16.2, 88.8;
    VectorXd curve1(2); curve1 << -18, 96;

    assert((curve.dEvaluate(2, 0) - curve0).norm() < 0.01);
    assert((curve.dEvaluate(2, 0.1) - curve01).norm() < 0.01);
    assert((curve.dEvaluate(2, 0.2) - curve02).norm() < 0.01);
    assert((curve.dEvaluate(2, 0.3) - curve03).norm() < 0.01);
    assert((curve.dEvaluate(2, 0.4) - curve04).norm() < 0.01);
    assert((curve.dEvaluate(2, 0.5) - curve05).norm() < 0.01);
    assert((curve.dEvaluate(2, 0.6) - curve06).norm() < 0.01);
    assert((curve.dEvaluate(2, 0.7) - curve07).norm() < 0.01);
    assert((curve.dEvaluate(2, 0.8) - curve08).norm() < 0.01);
    assert((curve.dEvaluate(2, 0.9) - curve09).norm() < 0.01);
    assert((curve.dEvaluate(2, 1) - curve1).norm() < 0.01);
}

void test_d_d_cubic_curve_evaluation_time_scaled() {    
    // initialize vectors
    VectorXd P0(2); P0 << 10, 5;
    VectorXd P1(2); P1 << 8, 2;
    VectorXd P2(2); P2 << 6, 3;
    VectorXd P3(2); P3 << 1, 20;

    Bezier::Curve<double> curve({P0, P1, P2, P3}, 10);

    // check output vectors, round to 3 decimal places; points evaluated using python script, see top of file
    VectorXd curve0(2); curve0 << 0, 24.0/100.0;
    VectorXd curve01(2); curve01 << -1.8/100.0, 31.2/100.0;
    VectorXd curve02(2); curve02 << -3.6/100.0, 38.4/100.0;
    VectorXd curve03(2); curve03 << -5.4/100.0, 45.6/100.0;
    VectorXd curve04(2); curve04 << -7.2/100.0, 52.8/100.0;
    VectorXd curve05(2); curve05 << -9.0/100.0, 60.0/100.0;
    VectorXd curve06(2); curve06 << -10.8/100.0, 67.2/100.0;
    VectorXd curve07(2); curve07 << -12.6/100.0, 74.4/100.0;
    VectorXd curve08(2); curve08 << -14.4/100.0, 81.6/100.0;
    VectorXd curve09(2); curve09 << -16.2/100.0, 88.8/100.0;
    VectorXd curve1(2); curve1 << -18.0/100.0, 96.0/100.0;

    double error = (curve.dEvaluate(2, 0) - curve0).norm();

    assert((curve.dEvaluate(2, 0) - curve0).norm() < 0.01);
    assert((curve.dEvaluate(2, 1) - curve01).norm() < 0.01);
    assert((curve.dEvaluate(2, 2) - curve02).norm() < 0.01);
    assert((curve.dEvaluate(2, 3) - curve03).norm() < 0.01);
    assert((curve.dEvaluate(2, 4) - curve04).norm() < 0.01);
    assert((curve.dEvaluate(2, 5) - curve05).norm() < 0.01);
    assert((curve.dEvaluate(2, 6) - curve06).norm() < 0.01);
    assert((curve.dEvaluate(2, 7) - curve07).norm() < 0.01);
    assert((curve.dEvaluate(2, 8) - curve08).norm() < 0.01);
    assert((curve.dEvaluate(2, 9) - curve09).norm() < 0.01);
    assert((curve.dEvaluate(2, 10) - curve1).norm() < 0.01);
}

void test_d_d_cubic_curve_evaluation_time_scaled_shifted() {    
    // initialize vectors
    VectorXd P0(2); P0 << 10, 5;
    VectorXd P1(2); P1 << 8, 2;
    VectorXd P2(2); P2 << 6, 3;
    VectorXd P3(2); P3 << 1, 20;

    Bezier::Curve<double> curve({P0, P1, P2, P3}, 10, 5);

    // check output vectors, round to 3 decimal places; points evaluated using python script, see top of file
    VectorXd curve0(2); curve0 << 0, 24.0/100.0;
    VectorXd curve01(2); curve01 << -1.8/100.0, 31.2/100.0;
    VectorXd curve02(2); curve02 << -3.6/100.0, 38.4/100.0;
    VectorXd curve03(2); curve03 << -5.4/100.0, 45.6/100.0;
    VectorXd curve04(2); curve04 << -7.2/100.0, 52.8/100.0;
    VectorXd curve05(2); curve05 << -9.0/100.0, 60.0/100.0;
    VectorXd curve06(2); curve06 << -10.8/100.0, 67.2/100.0;
    VectorXd curve07(2); curve07 << -12.6/100.0, 74.4/100.0;
    VectorXd curve08(2); curve08 << -14.4/100.0, 81.6/100.0;
    VectorXd curve09(2); curve09 << -16.2/100.0, 88.8/100.0;
    VectorXd curve1(2); curve1 << -18.0/100.0, 96.0/100.0;

    assert((curve.dEvaluate(2, 5) - curve0).norm() < 0.01);
    assert((curve.dEvaluate(2, 6) - curve01).norm() < 0.01);
    assert((curve.dEvaluate(2, 7) - curve02).norm() < 0.01);
    assert((curve.dEvaluate(2, 8) - curve03).norm() < 0.01);
    assert((curve.dEvaluate(2, 9) - curve04).norm() < 0.01);
    assert((curve.dEvaluate(2, 10) - curve05).norm() < 0.01);
    assert((curve.dEvaluate(2, 11) - curve06).norm() < 0.01);
    assert((curve.dEvaluate(2, 12) - curve07).norm() < 0.01);
    assert((curve.dEvaluate(2, 13) - curve08).norm() < 0.01);
    assert((curve.dEvaluate(2, 14) - curve09).norm() < 0.01);
    assert((curve.dEvaluate(2, 15) - curve1).norm() < 0.01);
}

int main(int argc, char const *argv[]) {
    test_d_cubic_curve_evaluation();
    test_d_d_cubic_curve_evaluation();
    test_d_d_cubic_curve_evaluation_time_scaled();
    test_d_d_cubic_curve_evaluation_time_scaled_shifted();
}
