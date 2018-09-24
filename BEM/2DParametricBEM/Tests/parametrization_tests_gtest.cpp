#include <cassert>
#include <cmath>
#include <stdlib.h>

#include <iostream>
#include <stdexcept>
#include <utility>

#include <Eigen/Dense>
#include "gtest/gtest.h"

#define _USE_MATH_DEFINES //for pi
#include "parametrized_line.hpp"
#include "parametrized_semi_circle.hpp"
#include "parametrized_fourier_sum.hpp"
#include "abstract_bem_space.hpp"
#include "discontinuous_space.hpp"
#include "continuous_space.hpp"

double eps = 1e-2; // A global threshold for error

//LineParametrizationTest is hieararcy name, Parametrization is a test in this hierarchy
TEST(LineParametrizationTest,Parametrization) {
  using Point = std::pair<double,double>;
  Eigen::Vector2d x1; x1 << 0,1; // Point (0,1)
  Eigen::Vector2d x2; x2 << 1,0; // Point (1,0)
  parametricbem2d::ParametrizedLine parametrization(x1,x2);
  // Test point (0.5,0.5) corresponding to t = 0.0
  Eigen::Vector2d testpoint;
  testpoint = parametrization(0.);
  EXPECT_NEAR(0.5,testpoint(0),eps);
  EXPECT_NEAR(0.5,testpoint(1),eps);
}

TEST(LineParametrizationTest,Derivative) {
  using Point = std::pair<double,double>;
  Eigen::Vector2d x1; x1 << 0,1; // Point (0,1)
  Eigen::Vector2d x2; x2 << 1,0; // Point (1,0)
  parametricbem2d::ParametrizedLine parametrization(x1,x2);
  // Test point (0.5,0.5) corresponding to t = 0.0
  Eigen::Vector2d testpoint;
  testpoint = parametrization(0.);
  // Derivative at the point corresponding to t = 0.3
  testpoint = parametrization.Derivative(0.3);
  EXPECT_NEAR(-1,testpoint(1)/testpoint(0),eps);
}

TEST(SemiCircleParametrizationTest,Parametrization) {
  // Parametrized semi-circle with unit radius
  parametricbem2d::ParametrizedSemiCircle parametrization;
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
  EXPECT_NEAR(M_PI,length,eps);
}

TEST(FourierSumParametrizationTest,Parametrization) {
  Eigen::MatrixXd a(2,1); //cosine coefficients ; N = 1
  Eigen::MatrixXd b(2,1); //sine coefficients ; N = 1
  a << 1.,
       0.;
  b << 0.,
       1.;
  // The fourier sum parametrization using these coefficients reduces to (cos(t),sin(t))
  parametricbem2d::ParametrizedFourierSum parametrization(a,b);
  double t = rand()/RAND_MAX;
  t = 2 * t - 1; //range -1 to 1
  Eigen::Vector2d randompoint = parametrization(t);
  EXPECT_NEAR(1,randompoint.norm(),eps);
}

TEST(InterfaceTest,ParameterRange) {
  double a,b;
  std::tie(a,b) = parametricbem2d::ParametrizedLine::ParameterRange();
  EXPECT_EQ(a,-1);
  EXPECT_EQ(b,1);
  std::tie(a,b) = parametricbem2d::ParametrizedFourierSum::ParameterRange();
  EXPECT_EQ(a,-1);
  EXPECT_EQ(b,1);
  std::tie(a,b) = parametricbem2d::ParametrizedSemiCircle::ParameterRange();
  EXPECT_EQ(a,-1);
  EXPECT_EQ(b,1);
}

TEST(InterfaceTest,IsWithinParameterRange) {
  EXPECT_EQ(true,parametricbem2d::ParametrizedLine::IsWithinParameterRange(0.99));
  EXPECT_EQ(false,parametricbem2d::ParametrizedFourierSum::IsWithinParameterRange(1.01));
  EXPECT_EQ(true,parametricbem2d::ParametrizedSemiCircle::IsWithinParameterRange(-0.99));
}

TEST(BemSpace,DiscontinuousSpace0) {
  parametricbem2d::AbstractBEMSpace *space = new parametricbem2d::DiscontinuousSpace<0>();
  using BasisFunctionPointer = parametricbem2d::AbstractBEMSpace::BasisFunctionPointer;
  int Q = space->getQ();
  EXPECT_EQ(Q,1);
  std::vector<BasisFunctionPointer> bases = space->getShapeFunctions();
  double t = rand()/RAND_MAX;
  t = 2 * t - 1; //range -1 to 1
  EXPECT_EQ(bases[0](t),1);
}

TEST(BemSpace,DiscontinuousSpace1) {
  parametricbem2d::AbstractBEMSpace *space = new parametricbem2d::DiscontinuousSpace<1>();
  using BasisFunctionPointer = parametricbem2d::AbstractBEMSpace::BasisFunctionPointer;
  int Q = space->getQ();
  EXPECT_EQ(Q,2);
  std::vector<BasisFunctionPointer> bases = space->getShapeFunctions();
  double t = rand()/RAND_MAX;
  t = 2 * t - 1; //range -1 to 1
  EXPECT_EQ(bases[0](t),0.5);
  EXPECT_EQ(bases[1](t),0.5*t);
}

TEST(BemSpace,ContinuousSpace0) {
  parametricbem2d::AbstractBEMSpace *space = new parametricbem2d::ContinuousSpace<0>();
  using BasisFunctionPointer = parametricbem2d::AbstractBEMSpace::BasisFunctionPointer;
  int Q = space->getQ();
  EXPECT_EQ(Q,1);
  std::vector<BasisFunctionPointer> bases = space->getShapeFunctions();
  double t = rand()/RAND_MAX;
  t = 2 * t - 1; //range -1 to 1
  EXPECT_EQ(bases[0](t),1);
}

TEST(BemSpace,ContinuousSpace1) {
  parametricbem2d::AbstractBEMSpace *space = new parametricbem2d::ContinuousSpace<1>();
  using BasisFunctionPointer = parametricbem2d::AbstractBEMSpace::BasisFunctionPointer;
  int Q = space->getQ();
  EXPECT_EQ(Q,2);
  std::vector<BasisFunctionPointer> bases = space->getShapeFunctions();
  double t = rand()/RAND_MAX;
  t = 2 * t - 1; //range -1 to 1
  EXPECT_EQ(bases[0](t),0.5*(1+t));
  EXPECT_EQ(bases[1](t),0.5*(1-t));
}

TEST(BemSpace,ContinuousSpace2) {
  const int p = 2;
  parametricbem2d::AbstractBEMSpace *space = new parametricbem2d::ContinuousSpace<p>();
  using BasisFunctionPointer = parametricbem2d::AbstractBEMSpace::BasisFunctionPointer;
  int Q = space->getQ();
  EXPECT_EQ(Q,3);
  std::vector<BasisFunctionPointer> bases = space->getShapeFunctions();
  double t = rand()/RAND_MAX;
  t = 2 * t - 1; //range -1 to 1
  EXPECT_EQ(bases[0](t),0.5*(1+t));
  EXPECT_EQ(bases[1](t),0.5*(1-t));
  EXPECT_EQ(bases[2](t),(1-t*t));
}

int main(int argc, char **argv) {
  srand(time(NULL));
  // run tests
  ::testing::InitGoogleTest(&argc,argv);
  return RUN_ALL_TESTS();
}
