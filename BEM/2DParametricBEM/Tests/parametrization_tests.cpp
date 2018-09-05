#include <cassert>
#include <cmath>
#include <stdlib.h>

#include <iostream>
#include <stdexcept>
#include <utility>

#include <Eigen/Dense>

#define _USE_MATH_DEFINES //for pi
#include "parametrized_line.hpp"
#include "parametrized_semi_circle.hpp"
#include "parametrized_fourier_sum.hpp"

double eps = 1e-2; // A global threshold for error

void test_parametrized_line() {
  using Point = std::pair<double,double>;
  Eigen::Vector2d x1; x1 << 0,1; // Point (0,1)
  Eigen::Vector2d x2; x2 << 1,0; // Point (1,0)
  ParametrizedLine parametrization(x1,x2);
  // Test point (0.5,0.5) corresponding to t = 0.0
  Eigen::Vector2d testpoint;
  testpoint = parametrization(0.);
  assert(fabs(testpoint(0)-0.5)<eps);
  assert(fabs(testpoint(1)-0.5)<eps);
  std::cout << "(" << testpoint(0) << "," << testpoint(1) <<") ; expected (0.5,0.5)" << std::endl;
  // Derivative at the point corresponding to t = 0.3
  testpoint = parametrization.Derivative(0.3);
  assert(fabs(testpoint(1)/testpoint(0)+1)<eps);
  std::cout << "slope : "<<testpoint(1)/testpoint(0) <<" ; expected : -1" << std::endl;
}

void test_parametrized_semi_circle() {
  // Parametrized semi-circle with unit radius
  ParametrizedSemiCircle parametrization;
  double tmin,tmax;
  std::tie(tmin,tmax) = parametrization.ParameterRange();
  double length = 0.;
  double t1,t2;
  // Representing the semi-circular parametrized curve by Nbins number of
  // Points or Nbins -1 number of line segments and summing up the lengths
  // to get an approximation of Pi
  int Nbins = 1000;
  for (int i = 0 ; i<Nbins-1 ; ++i) {
    t1 = tmin + i * (tmax-tmin)/(Nbins-1); // Current Point
    t2 = tmin + (i+1) * (tmax-tmin)/(Nbins-1); // Next Point
    Eigen::Vector2d start,end;
    start = parametrization(t1);
    end = parametrization(t2);
    length += (end-start).norm(); // Adding the length of current line segment
  }
  std::cout << "Approximate value of Pi: " << length << std::endl;
  assert(fabs(length-M_PI)<eps);
}

void test_fourier_parametrization() {
  Eigen::MatrixXd a(2,1); //cosine coefficients ; N = 1
  Eigen::MatrixXd b(2,1); //sine coefficients ; N = 1
  a << 1.,
       0.;
  b << 0.,
       1.;
  // The fourier sum parametrization using these coefficients reduces to (cos(t),sin(t))
  ParametrizedFourierSum parametrization(a,b);
  double t = rand()/RAND_MAX;
  t = 2 * t - 1; //range -1 to 1
  Eigen::Vector2d randompoint = parametrization(t);
  std::cout << "Norm at random point = " << randompoint.norm() << " (expected value = 1.)" << std::endl;
  assert(fabs(randompoint.norm()-1)<eps);
}

int main() {
  srand(time(NULL));
  // run tests
  test_parametrized_line();
  test_parametrized_semi_circle();
  test_fourier_parametrization();
  return 0;
}
