#include <cassert>
#include <cmath>
#include <iostream>
#include <stdexcept>
#include <stdlib.h>
#include <utility>

#include "BoundaryMesh.hpp"
#include "abstract_bem_space.hpp"
#include "buildK.hpp"
#include "buildV.hpp"
#include "continuous_space.hpp"
#include "discontinuous_space.hpp"
#include "doubleLayerPotential.hpp"
#include "double_layer.hpp"
#include "integral_gauss.hpp"
#include "parametrized_circular_arc.hpp"
#include "parametrized_fourier_sum.hpp"
#include "parametrized_line.hpp"
#include "parametrized_mesh.hpp"
#include "parametrized_semi_circle.hpp"
#include "singleLayerPotential.hpp"
#include "single_layer.hpp"
#include "gtest/gtest.h"
#include <Eigen/Dense>

#define _USE_MATH_DEFINES // for pi

double eps = 1e-5; // A global threshold for error

// LineParametrizationTest is hieararcy name, Parametrization is a test in this
// hierarchy
TEST(LineParametrizationTest, Parametrization) {
  using Point = std::pair<double, double>;
  Eigen::Vector2d x1;
  x1 << 0, 1; // Point (0,1)
  Eigen::Vector2d x2;
  x2 << 1, 0; // Point (1,0)
  parametricbem2d::ParametrizedLine parametrization(x1, x2);
  // Test point (0.5,0.5) corresponding to t = 0.0
  Eigen::Vector2d testpoint;
  testpoint = parametrization(0.);
  EXPECT_NEAR(0.5, testpoint(0), eps);
  EXPECT_NEAR(0.5, testpoint(1), eps);
}

TEST(LineParametrizationTest, Derivative) {
  using Point = std::pair<double, double>;
  Eigen::Vector2d x1;
  x1 << 0, 1; // Point (0,1)
  Eigen::Vector2d x2;
  x2 << 1, 0; // Point (1,0)
  parametricbem2d::ParametrizedLine parametrization(x1, x2);
  // Test point (0.5,0.5) corresponding to t = 0.0
  Eigen::Vector2d testpoint;
  testpoint = parametrization(0.);
  // Derivative at the point corresponding to t = 0.3
  testpoint = parametrization.Derivative(0.3);
  EXPECT_NEAR(-1, testpoint(1) / testpoint(0), eps);
}

TEST(SemiCircleParametrizationTest, Parametrization) {
  // Parametrized semi-circle with unit radius
  parametricbem2d::ParametrizedSemiCircle parametrization;
  double tmin, tmax;
  std::tie(tmin, tmax) = parametrization.ParameterRange();
  double length = 0.;
  double t1, t2;
  // Representing the semi-circular parametrized curve by Nbins number of
  // Points or Nbins -1 number of line segments and summing up the lengths
  // to get an approximation of Pi
  int Nbins = 1000;
  for (int i = 0; i < Nbins - 1; ++i) {
    t1 = tmin + i * (tmax - tmin) / (Nbins - 1);       // Current Point
    t2 = tmin + (i + 1) * (tmax - tmin) / (Nbins - 1); // Next Point
    Eigen::Vector2d start, end;
    start = parametrization(t1);
    end = parametrization(t2);
    length += (end - start).norm(); // Adding the length of current line segment
  }
  EXPECT_NEAR(M_PI, length, eps);
}

TEST(ParametrizedCircularArcTest, Parametrization) {
  // Parametrized with unit radius, centered at (1,-1)
  Eigen::Vector2d center(2);
  center << 1, -1;
  double r = 1.2;
  parametricbem2d::ParametrizedCircularArc parametrization(center, r, 0,
                                                           M_PI / 2);
  double tmin, tmax;
  std::tie(tmin, tmax) = parametrization.ParameterRange();
  double length = 0.;
  double t1, t2;
  // Representing the circular arc parametrized curve by Nbins number of
  // Points or Nbins -1 number of line segments and summing up the lengths
  // to get an approximation of Pi
  int Nbins = 1000;
  for (int i = 0; i < Nbins - 1; ++i) {
    t1 = tmin + i * (tmax - tmin) / (Nbins - 1);       // Current Point
    t2 = tmin + (i + 1) * (tmax - tmin) / (Nbins - 1); // Next Point
    Eigen::Vector2d start, end;
    start = parametrization(t1);
    end = parametrization(t2);
    length += (end - start).norm(); // Adding the length of current line segment
  }
  EXPECT_NEAR(M_PI, length * 2 / r, eps);
}

TEST(FourierSumParametrizationTest, Parametrization) {
  Eigen::MatrixXd a(2, 1); // cosine coefficients ; N = 1
  Eigen::MatrixXd b(2, 1); // sine coefficients ; N = 1
  a << 1., 0.;
  b << 0., 1.;
  // The fourier sum parametrization using these coefficients reduces to
  // (cos(t),sin(t))
  parametricbem2d::ParametrizedFourierSum parametrization(a, b);
  double t = rand() / RAND_MAX;
  t = 2 * t - 1; // range -1 to 1
  Eigen::Vector2d randompoint = parametrization(t);
  EXPECT_NEAR(1, randompoint.norm(), eps);
}

TEST(InterfaceTest, ParameterRange) {
  double a, b;
  std::tie(a, b) = parametricbem2d::ParametrizedLine::ParameterRange();
  EXPECT_EQ(a, -1);
  EXPECT_EQ(b, 1);
  std::tie(a, b) = parametricbem2d::ParametrizedFourierSum::ParameterRange();
  EXPECT_EQ(a, -1);
  EXPECT_EQ(b, 1);
  std::tie(a, b) = parametricbem2d::ParametrizedSemiCircle::ParameterRange();
  EXPECT_EQ(a, -1);
  EXPECT_EQ(b, 1);
}

TEST(InterfaceTest, IsWithinParameterRange) {
  EXPECT_EQ(true,
            parametricbem2d::ParametrizedLine::IsWithinParameterRange(0.99));
  EXPECT_EQ(
      false,
      parametricbem2d::ParametrizedFourierSum::IsWithinParameterRange(1.01));
  EXPECT_EQ(
      true,
      parametricbem2d::ParametrizedSemiCircle::IsWithinParameterRange(-0.99));
}

double integrand(double x) { return x * x; }

TEST(IntegralGauss, ComputeIntegral) {
  double integral = parametricbem2d::ComputeIntegral(integrand, 0, 1, 10);
  EXPECT_NEAR(integral, 1. / 3., eps);
}

TEST(BemSpace, DiscontinuousSpace0) {
  parametricbem2d::AbstractBEMSpace *space =
      new parametricbem2d::DiscontinuousSpace<0>();
  int Q = space->getQ();
  EXPECT_EQ(Q, 1);
  double t = rand() / RAND_MAX;
  t = 2 * t - 1; // range -1 to 1
  EXPECT_EQ(space->evaluateShapeFunction(0,t), 1);
  delete space;
}

TEST(BemSpace, DiscontinuousSpace1) {
  parametricbem2d::AbstractBEMSpace *space =
      new parametricbem2d::DiscontinuousSpace<1>();
  int Q = space->getQ();
  EXPECT_EQ(Q, 2);
  double t = rand() / RAND_MAX;
  t = 2 * t - 1; // range -1 to 1
  EXPECT_EQ(space->evaluateShapeFunction(0,t), 0.5);
  EXPECT_EQ(space->evaluateShapeFunction(1,t), 0.5 * t);
  delete space;
}

TEST(BemSpace, ContinuousSpace1) {
  parametricbem2d::AbstractBEMSpace *space =
      new parametricbem2d::ContinuousSpace<1>();
  int Q = space->getQ();
  EXPECT_EQ(Q, 2);
  double t = rand() / RAND_MAX;
  t = 2 * t - 1; // range -1 to 1
  EXPECT_EQ(space->evaluateShapeFunction(0,t), 0.5 * (1 + t));
  EXPECT_EQ(space->evaluateShapeFunction(1,t), 0.5 * (1 - t));
  delete space;
}

TEST(BemSpace, ContinuousSpace2) {
  const int p = 2;
  parametricbem2d::AbstractBEMSpace *space =
      new parametricbem2d::ContinuousSpace<p>();
  int Q = space->getQ();
  EXPECT_EQ(Q, 3);
  double t = rand() / RAND_MAX;
  t = 2 * t - 1; // range -1 to 1
  EXPECT_EQ(space->evaluateShapeFunction(0,t), 0.5 * (1 + t));
  EXPECT_EQ(space->evaluateShapeFunction(1,t), 0.5 * (1 - t));
  EXPECT_EQ(space->evaluateShapeFunction(2,t), (1 - t * t));
  delete space;
}

TEST(LocGlobMap, ContinuousSpace1) {
  parametricbem2d::AbstractBEMSpace *space =
      new parametricbem2d::ContinuousSpace<1>();
  EXPECT_EQ(space->LocGlobMap(1, 1, 4), 2);
  EXPECT_EQ(space->LocGlobMap(1, 2, 4), 3);
  EXPECT_EQ(space->LocGlobMap(1, 3, 4), 4);
  EXPECT_EQ(space->LocGlobMap(2, 3, 4), 3);
  delete space;
}

TEST(LocGlobMap, DiscontinuousSpace0) {
  parametricbem2d::AbstractBEMSpace *space =
      new parametricbem2d::DiscontinuousSpace<0>();
  EXPECT_EQ(space->LocGlobMap(1, 1, 4), 1);
  EXPECT_EQ(space->LocGlobMap(1, 2, 4), 2);
  delete space;
}

TEST(LocGlobMap, DiscontinuousSpace1) {
  parametricbem2d::AbstractBEMSpace *space =
      new parametricbem2d::DiscontinuousSpace<1>();
  EXPECT_EQ(space->LocGlobMap(2, 1, 4), 5);
  EXPECT_EQ(space->LocGlobMap(2, 2, 4), 6);
  EXPECT_EQ(space->LocGlobMap(1, 4, 4), 4);
  EXPECT_EQ(space->LocGlobMap(2, 4, 4), 8);
  delete space;
}

TEST(Split, ParametrizedLine) {
  using Point = std::pair<double, double>;
  using PanelVector = parametricbem2d::PanelVector;
  Eigen::Vector2d x1;
  x1 << 0, 1; // Point (0,1)
  Eigen::Vector2d x2;
  x2 << 1, 0; // Point (1,0)
  parametricbem2d::ParametrizedLine parametrization(x1, x2);
  unsigned int N = 10;
  PanelVector components = parametrization.split(N);
  for (unsigned int i = 0; i < N - 1; ++i) {
    EXPECT_NEAR(
        (components[i]->operator()(1) - components[i + 1]->operator()(-1))
            .norm(),
        0, eps);
  }
}

TEST(Split, ParametrizedCircularArc) {
  using Point = std::pair<double, double>;
  using PanelVector = parametricbem2d::PanelVector;
  Eigen::Vector2d center;
  center << 2.5, 3.6;
  double radius = 33;
  parametricbem2d::ParametrizedCircularArc parametrization(
      center, radius, .33 * M_PI, .99 * M_PI);
  unsigned int N = 10;
  PanelVector components = parametrization.split(N);
  for (unsigned int i = 0; i < N - 1; ++i) {
    // std::cout << "i: " << i << std::endl;
    EXPECT_NEAR(
        (components[i]->operator()(1) - components[i + 1]->operator()(-1))
            .norm(),
        0, eps);
  }
}

TEST(Split, ParametrizedFourierSum) {
  Eigen::MatrixXd a(2, 1); // cosine coefficients ; N = 1
  Eigen::MatrixXd b(2, 1); // sine coefficients ; N = 1
  a << 1., 0.;
  b << 0., 1.;
  parametricbem2d::ParametrizedFourierSum parametrization(a, b);
  unsigned int N = 10;
  using PanelVector = parametricbem2d::PanelVector;
  PanelVector components = parametrization.split(N);
  for (unsigned int i = 0; i < N - 1; ++i)
    EXPECT_NEAR(
        (components[i]->operator()(1) - components[i + 1]->operator()(-1))
            .norm(),
        0, eps);
}

TEST(ParametrizedMeshTest, MemberFunctions) {
  using PanelVector = parametricbem2d::PanelVector;
  Eigen::RowVectorXd x1(2);
  x1 << 0, 0; // Point (0,0)
  Eigen::RowVectorXd x2(2);
  x2 << 1, 0; // Point (1,0)
  Eigen::RowVectorXd x3(2);
  x3 << 1, .5; // Point (1,0.5)
  Eigen::RowVectorXd x4(2);
  x4 << 0, 1.5; // Point (0,1.5)
  parametricbem2d::ParametrizedLine line1(x1, x2);
  parametricbem2d::ParametrizedLine line2(x2, x3);
  parametricbem2d::ParametrizedLine line3(x3, x4);
  parametricbem2d::ParametrizedLine line4(x4, x1);
  PanelVector line1panels = line1.split(2);
  PanelVector line2panels = line2.split(2);
  PanelVector line3panels = line3.split(2);
  PanelVector line4panels = line4.split(2);
  PanelVector panels;
  panels.insert(panels.end(), line1panels.begin(), line1panels.end());
  panels.insert(panels.end(), line2panels.begin(), line2panels.end());
  panels.insert(panels.end(), line3panels.begin(), line3panels.end());
  panels.insert(panels.end(), line4panels.begin(), line4panels.end());
  parametricbem2d::ParametrizedMesh mesh(panels);
  Eigen::Vector2d vertex2 = mesh.getVertex(1);
  EXPECT_EQ(mesh.getNumPanels(), 8);
  EXPECT_NEAR(vertex2(0), 0.5, eps);
  EXPECT_NEAR(vertex2(1), 0, eps);
}

TEST(SingleLayer_0, PanelOrientedAssembly) {
  using PanelVector = parametricbem2d::PanelVector;
  Eigen::RowVectorXd x1(2);
  x1 << 0, 0; // Point (0,0)
  Eigen::RowVectorXd x2(2);
  x2 << 1, 0; // Point (1,0)
  Eigen::RowVectorXd x3(2);
  x3 << 1, .5; // Point (1,0.5)
  Eigen::RowVectorXd x4(2);
  x4 << 0, 1.5; // Point (0,1.5)
  parametricbem2d::ParametrizedLine line1(x1, x2);
  parametricbem2d::ParametrizedLine line2(x2, x3);
  parametricbem2d::ParametrizedLine line3(x3, x4);
  parametricbem2d::ParametrizedLine line4(x4, x1);
  PanelVector line1panels = line1.split(1);
  PanelVector line2panels = line2.split(2);
  PanelVector line3panels = line3.split(3);
  PanelVector line4panels = line4.split(4);
  PanelVector panels;
  panels.insert(panels.end(), line1panels.begin(), line1panels.end());
  panels.insert(panels.end(), line2panels.begin(), line2panels.end());
  panels.insert(panels.end(), line3panels.begin(), line3panels.end());
  panels.insert(panels.end(), line4panels.begin(), line4panels.end());
  parametricbem2d::ParametrizedMesh mesh(panels);
  parametricbem2d::DiscontinuousSpace<0> space;
  Eigen::MatrixXd galerkin =
      parametricbem2d::single_layer::SingleLayerMatrix(mesh, space, 32);
  unsigned int numpanels = mesh.getNumPanels();
  EXPECT_EQ(numpanels, galerkin.cols());
  EXPECT_EQ(numpanels, galerkin.rows());
}

TEST(SingleLayer_1, PanelOrientedAssembly) {
  using PanelVector = parametricbem2d::PanelVector;
  Eigen::RowVectorXd x1(2);
  x1 << 0, 0; // Point (0,0)
  Eigen::RowVectorXd x2(2);
  x2 << 1, 0; // Point (1,0)
  Eigen::RowVectorXd x3(2);
  x3 << 1, .5; // Point (1,0.5)
  Eigen::RowVectorXd x4(2);
  x4 << 0, 1.5; // Point (0,1.5)
  parametricbem2d::ParametrizedLine line1(x1, x2);
  parametricbem2d::ParametrizedLine line2(x2, x3);
  parametricbem2d::ParametrizedLine line3(x3, x4);
  parametricbem2d::ParametrizedLine line4(x4, x1);
  PanelVector line1panels = line1.split(1);
  PanelVector line2panels = line2.split(1);
  PanelVector line3panels = line3.split(1);
  PanelVector line4panels = line4.split(1);
  PanelVector panels;
  panels.insert(panels.end(), line1panels.begin(), line1panels.end());
  panels.insert(panels.end(), line2panels.begin(), line2panels.end());
  panels.insert(panels.end(), line3panels.begin(), line3panels.end());
  panels.insert(panels.end(), line4panels.begin(), line4panels.end());
  parametricbem2d::ParametrizedMesh mesh(panels);
  parametricbem2d::DiscontinuousSpace<1> space;
  Eigen::MatrixXd galerkin =
      parametricbem2d::single_layer::SingleLayerMatrix(mesh, space, 32);
  unsigned int numpanels = mesh.getNumPanels();
  EXPECT_EQ(8, galerkin.cols());
  EXPECT_EQ(8, galerkin.rows());
}

TEST(SingleLayer_0, CppHilbertComparison) {
  using PanelVector = parametricbem2d::PanelVector;
  Eigen::RowVectorXd x1(2);
  x1 << 0, 0; // Point (0,0)
  Eigen::RowVectorXd x2(2);
  x2 << 1, 0; // Point (1,0)
  Eigen::RowVectorXd x3(2);
  x3 << 1, .5; // Point (1,0.5)
  Eigen::RowVectorXd x4(2);
  x4 << 0, 1.5; // Point (0,1.5)
  parametricbem2d::ParametrizedLine line1(x1, x2);
  parametricbem2d::ParametrizedLine line2(x2, x3);
  parametricbem2d::ParametrizedLine line3(x3, x4);
  parametricbem2d::ParametrizedLine line4(x4, x1);
  PanelVector line1panels = line1.split(1);
  PanelVector line2panels = line2.split(1);
  PanelVector line3panels = line3.split(1);
  PanelVector line4panels = line4.split(1);
  PanelVector panels;
  panels.insert(panels.end(), line1panels.begin(), line1panels.end());
  panels.insert(panels.end(), line2panels.begin(), line2panels.end());
  panels.insert(panels.end(), line3panels.begin(), line3panels.end());
  panels.insert(panels.end(), line4panels.begin(), line4panels.end());
  parametricbem2d::ParametrizedMesh parametrizedmesh(panels);
  parametricbem2d::DiscontinuousSpace<0> space;
  Eigen::MatrixXd galerkinnew =
      parametricbem2d::single_layer::SingleLayerMatrix(parametrizedmesh, space,
                                                       32);
  Eigen::MatrixXd coords(4, 2);
  coords << x1, x2, x3, x4;
  Eigen::Matrix<int, 4, 2> elems;
  elems << 0, 1, 1, 2, 2, 3, 3, 0;
  BoundaryMesh boundarymesh(coords, elems);
  Eigen::MatrixXd galerkinold;
  computeV(galerkinold, boundarymesh, 0);
  EXPECT_NEAR((galerkinold - galerkinnew).norm(), 0, eps);
}

TEST(DoubleLayer_1_0, PanelOrientedAssembly) {
  using PanelVector = parametricbem2d::PanelVector;
  Eigen::RowVectorXd x1(2);
  x1 << 0, 0; // Point (0,0)
  Eigen::RowVectorXd x2(2);
  x2 << 1, 0; // Point (1,0)
  Eigen::RowVectorXd x3(2);
  x3 << 1, .5; // Point (1,0.5)
  Eigen::RowVectorXd x4(2);
  x4 << 0, 1.5; // Point (0,1.5)
  parametricbem2d::ParametrizedLine line1(x1, x2);
  parametricbem2d::ParametrizedLine line2(x2, x3);
  parametricbem2d::ParametrizedLine line3(x3, x4);
  parametricbem2d::ParametrizedLine line4(x4, x1);
  PanelVector line1panels = line1.split(1);
  PanelVector line2panels = line2.split(2);
  PanelVector line3panels = line3.split(3);
  PanelVector line4panels = line4.split(4);
  PanelVector panels;
  panels.insert(panels.end(), line1panels.begin(), line1panels.end());
  panels.insert(panels.end(), line2panels.begin(), line2panels.end());
  panels.insert(panels.end(), line3panels.begin(), line3panels.end());
  panels.insert(panels.end(), line4panels.begin(), line4panels.end());
  parametricbem2d::ParametrizedMesh mesh(panels);
  parametricbem2d::DiscontinuousSpace<0> test_space;
  parametricbem2d::ContinuousSpace<1> trial_space;
  Eigen::MatrixXd galerkin = parametricbem2d::double_layer::DoubleLayerMatrix(
      mesh, trial_space, test_space, 32);
  unsigned int numpanels = mesh.getNumPanels();
  EXPECT_EQ(numpanels, galerkin.cols());
  EXPECT_EQ(numpanels, galerkin.rows());
}

TEST(DoubleLayer_1_0, CppHilbertComparison) {
  using PanelVector = parametricbem2d::PanelVector;
  Eigen::RowVectorXd x1(2);
  x1 << 0, 0; // Point (0,0)
  Eigen::RowVectorXd x2(2);
  x2 << 1, 0; // Point (1,0)
  Eigen::RowVectorXd x3(2);
  x3 << 1, .5; // Point (1,0.5)
  Eigen::RowVectorXd x4(2);
  x4 << 0, 1.5; // Point (0,1.5)
  parametricbem2d::ParametrizedLine line1(x1, x2);
  parametricbem2d::ParametrizedLine line2(x2, x3);
  parametricbem2d::ParametrizedLine line3(x3, x4);
  parametricbem2d::ParametrizedLine line4(x4, x1);
  PanelVector line1panels = line1.split(1);
  PanelVector line2panels = line2.split(1);
  PanelVector line3panels = line3.split(1);
  PanelVector line4panels = line4.split(1);
  PanelVector panels;
  panels.insert(panels.end(), line1panels.begin(), line1panels.end());
  panels.insert(panels.end(), line2panels.begin(), line2panels.end());
  panels.insert(panels.end(), line3panels.begin(), line3panels.end());
  panels.insert(panels.end(), line4panels.begin(), line4panels.end());
  parametricbem2d::ParametrizedMesh parametrizedmesh(panels);
  parametricbem2d::DiscontinuousSpace<0> test_space;
  parametricbem2d::ContinuousSpace<1> trial_space;
  Eigen::MatrixXd galerkin = parametricbem2d::double_layer::DoubleLayerMatrix(
      parametrizedmesh, trial_space, test_space, 32);
  Eigen::MatrixXd coords(4, 2);
  coords << x1, x2, x3, x4;
  Eigen::Matrix<int, 4, 2> elems;
  elems << 0, 1, 1, 2, 2, 3, 3, 0;
  BoundaryMesh boundarymesh(coords, elems);
  Eigen::MatrixXd galerkinold;
  computeK(galerkinold, boundarymesh, 0);
  EXPECT_NEAR((galerkinold - galerkin).norm(), 0, eps);
}

TEST(DoubleLayer_0_0, PanelOrientedAssembly) {
  using PanelVector = parametricbem2d::PanelVector;
  Eigen::RowVectorXd x1(2);
  x1 << 0, 0; // Point (0,0)
  Eigen::RowVectorXd x2(2);
  x2 << 1, 0; // Point (1,0)
  Eigen::RowVectorXd x3(2);
  x3 << 1, .5; // Point (1,0.5)
  Eigen::RowVectorXd x4(2);
  x4 << 0, 1.5; // Point (0,1.5)
  parametricbem2d::ParametrizedLine line1(x1, x2);
  parametricbem2d::ParametrizedLine line2(x2, x3);
  parametricbem2d::ParametrizedLine line3(x3, x4);
  parametricbem2d::ParametrizedLine line4(x4, x1);
  PanelVector line1panels = line1.split(1);
  PanelVector line2panels = line2.split(2);
  PanelVector line3panels = line3.split(3);
  PanelVector line4panels = line4.split(4);
  PanelVector panels;
  panels.insert(panels.end(), line1panels.begin(), line1panels.end());
  panels.insert(panels.end(), line2panels.begin(), line2panels.end());
  panels.insert(panels.end(), line3panels.begin(), line3panels.end());
  panels.insert(panels.end(), line4panels.begin(), line4panels.end());
  parametricbem2d::ParametrizedMesh mesh(panels);
  parametricbem2d::DiscontinuousSpace<0> test_space;
  parametricbem2d::DiscontinuousSpace<0> trial_space;
  Eigen::MatrixXd galerkin = parametricbem2d::double_layer::DoubleLayerMatrix(
      mesh, trial_space, test_space, 32);
  unsigned int numpanels = mesh.getNumPanels();
  EXPECT_EQ(numpanels, galerkin.cols());
  EXPECT_EQ(numpanels, galerkin.rows());
}

TEST(DoubleLayer_0_0, CppHilbertComparison) {
  using PanelVector = parametricbem2d::PanelVector;
  Eigen::RowVectorXd x1(2);
  x1 << 0, 0; // Point (0,0)
  Eigen::RowVectorXd x2(2);
  x2 << 1, 0; // Point (1,0)
  Eigen::RowVectorXd x3(2);
  x3 << 1, .5; // Point (1,0.5)
  Eigen::RowVectorXd x4(2);
  x4 << 0, 1.5; // Point (0,1.5)
  parametricbem2d::ParametrizedLine line1(x1, x2);
  parametricbem2d::ParametrizedLine line2(x2, x3);
  parametricbem2d::ParametrizedLine line3(x3, x4);
  parametricbem2d::ParametrizedLine line4(x4, x1);
  PanelVector line1panels = line1.split(1);
  PanelVector line2panels = line2.split(1);
  PanelVector line3panels = line3.split(1);
  PanelVector line4panels = line4.split(1);
  PanelVector panels;
  panels.insert(panels.end(), line1panels.begin(), line1panels.end());
  panels.insert(panels.end(), line2panels.begin(), line2panels.end());
  panels.insert(panels.end(), line3panels.begin(), line3panels.end());
  panels.insert(panels.end(), line4panels.begin(), line4panels.end());
  parametricbem2d::ParametrizedMesh parametrizedmesh(panels);
  parametricbem2d::DiscontinuousSpace<0> test_space;
  parametricbem2d::DiscontinuousSpace<0> trial_space;
  Eigen::MatrixXd galerkin = parametricbem2d::double_layer::DoubleLayerMatrix(
      parametrizedmesh, trial_space, test_space, 32);
  Eigen::MatrixXd coords(4, 2);
  coords << x1, x2, x3, x4;
  Eigen::Matrix<int, 4, 2> elems;
  elems << 0, 1, 1, 2, 2, 3, 3, 0;
  BoundaryMesh boundarymesh(coords, elems);
  Eigen::MatrixXd galerkinold;
  computeK00(galerkinold, boundarymesh, 0);
  EXPECT_NEAR((galerkinold - galerkin).norm(), 0, eps);
}

int main(int argc, char **argv) {
  srand(time(NULL));
  // run tests
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
