#include <cassert>
#include <cmath>
#include <stdlib.h>
#include <iostream>
#include <stdexcept>
#include <utility>

#include <Eigen/Dense>
#include "gtest/gtest.h"
#include "parametrized_line.hpp"
#include "parametrized_semi_circle.hpp"
#include "parametrized_fourier_sum.hpp"
#include "parametrized_circular_arc.hpp"
#include "abstract_bem_space.hpp"
#include "discontinuous_space.hpp"
#include "continuous_space.hpp"
#include "single_layer.hpp"
#include "singleLayerPotential.hpp"
#include "parametrized_mesh.hpp"
#include "BoundaryMesh.hpp"
#include "buildV.hpp"


#define _USE_MATH_DEFINES //for pi

double eps = 1e-5; // A global threshold for error

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

TEST(ParametrizedCircularArcTest,Parametrization) {
  // Parametrized with unit radius, centered at (1,-1)
  Eigen::Vector2d center(2); center << 1,-1;
  double r = 1.2;
  parametricbem2d::ParametrizedCircularArc parametrization(center,r,0,M_PI/2);
  double tmin,tmax;
  std::tie(tmin,tmax) = parametrization.ParameterRange();
  double length = 0.;
  double t1,t2;
  // Representing the circular arc parametrized curve by Nbins number of
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
  EXPECT_NEAR(M_PI,length*2/r,eps);
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
  delete space;
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
  delete space;
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
  delete space;
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
  delete space;
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
  delete space;
}

TEST(SingleLayer_0,CoincidingPanels) {
  // Test for the single layer
  Eigen::Vector2d x1; x1 << 10,9; // Point (0,1)
  Eigen::Vector2d x2; x2 << 7,11; // Point (1,0)
  parametricbem2d::ParametrizedLine parametrization(x1,x2);
  //parametricbem2d::AbstractBEMSpace *space = new parametricbem2d::DiscontinuousSpace<0>();
  parametricbem2d::DiscontinuousSpace<0> space;
  Eigen::MatrixXd interaction_matrix = SingleLayer(parametrization,
                                                   parametrization,
                                                   space,
                                                   32);

  double old_lib_soln = computeVij(x1,x2,x1,x2,0.);
  EXPECT_NEAR(interaction_matrix(0,0),old_lib_soln,eps);
}

TEST(SingleLayer_0,AdjacentPanels) {
  // Test for the single layer
  Eigen::Vector2d x1; x1 << 3,4; // Point (0,1)
  Eigen::Vector2d x2; x2 << 7,8; // Point (1,0)
  Eigen::Vector2d x3; x3 << 11,15; // Point (0,-1)
  parametricbem2d::ParametrizedLine panel1(x1,x2);
  parametricbem2d::ParametrizedLine panel2(x2,x3);
  //parametricbem2d::AbstractBEMSpace *space = new parametricbem2d::DiscontinuousSpace<0>();
  parametricbem2d::DiscontinuousSpace<0> space;
  Eigen::MatrixXd interaction_matrix = SingleLayer(panel1,
                                                   panel2,
                                                   space,
                                                   32);

  double old_lib_soln = computeVij(x1,x2,x2,x3,0.);
  EXPECT_NEAR(interaction_matrix(0,0),old_lib_soln,eps);
}

TEST(SingleLayer_0,General) {
  // Test for the single layer
  Eigen::Vector2d x1; x1 << 0,1; // Point (0,1)
  Eigen::Vector2d x2; x2 << 1,0; // Point (1,0)
  Eigen::Vector2d x3; x3 << 0,-1; // Point (0,-1)
  Eigen::Vector2d x4; x4 << -1,0; // Point (-1,0)
  parametricbem2d::ParametrizedLine panel1(x1,x2);
  parametricbem2d::ParametrizedLine panel2(x3,x4);
  parametricbem2d::DiscontinuousSpace<0> space;
  Eigen::MatrixXd interaction_matrix = SingleLayer(panel1,
                                                   panel2,
                                                   space,
                                                   33);

  double old_lib_soln = computeVij(x1,x2,x3,x4,0.);
  EXPECT_NEAR(interaction_matrix(0,0),old_lib_soln,eps);
}

TEST(LocGlobMap,ContinuousSpace0) {
  parametricbem2d::AbstractBEMSpace *space = new parametricbem2d::ContinuousSpace<0>();
  using LocGlobMapPointer = parametricbem2d::AbstractBEMSpace::LocGlobMapPointer;
  LocGlobMapPointer map = space->getLocGlobMap();
  EXPECT_EQ(map(1,5,10),1);
  delete space;
}

TEST(LocGlobMap,ContinuousSpace1) {
  parametricbem2d::AbstractBEMSpace *space = new parametricbem2d::ContinuousSpace<1>();
  using LocGlobMapPointer = parametricbem2d::AbstractBEMSpace::LocGlobMapPointer;
  LocGlobMapPointer map = space->getLocGlobMap();
  EXPECT_EQ(map(1,1,4),1);
  EXPECT_EQ(map(1,2,4),2);
  EXPECT_EQ(map(1,1,4),1);
  EXPECT_EQ(map(2,1,4),4);
  delete space;
}

TEST(LocGlobMap,DiscontinuousSpace0) {
  parametricbem2d::AbstractBEMSpace *space = new parametricbem2d::DiscontinuousSpace<0>();
  using LocGlobMapPointer = parametricbem2d::AbstractBEMSpace::LocGlobMapPointer;
  LocGlobMapPointer map = space->getLocGlobMap();
  EXPECT_EQ(map(1,1,4),1);
  EXPECT_EQ(map(1,2,4),2);
  delete space;
}

TEST(LocGlobMap,DiscontinuousSpace1) {
  parametricbem2d::AbstractBEMSpace *space = new parametricbem2d::DiscontinuousSpace<1>();
  using LocGlobMapPointer = parametricbem2d::AbstractBEMSpace::LocGlobMapPointer;
  LocGlobMapPointer map = space->getLocGlobMap();
  EXPECT_EQ(map(2,1,4),5);
  EXPECT_EQ(map(2,2,4),6);
  EXPECT_EQ(map(1,4,4),4);
  EXPECT_EQ(map(2,4,4),8);
  delete space;
}

TEST(Split,ParametrizedLine) {
  using Point = std::pair<double,double>;
  using PanelVector = parametricbem2d::PanelVector;
  Eigen::Vector2d x1; x1 << 0,1; // Point (0,1)
  Eigen::Vector2d x2; x2 << 1,0; // Point (1,0)
  parametricbem2d::ParametrizedLine parametrization(x1,x2);
  unsigned int N = 10;
  PanelVector components = parametrization.split(N);
  for (unsigned int i = 0 ; i < N-1 ; ++i) {
    EXPECT_NEAR((components[i]->operator()(1) -
                 components[i+1]->operator()(-1)).norm(),0,eps);
  }
}

TEST(Split,ParametrizedCircularArc) {
  using Point = std::pair<double,double>;
  using PanelVector = parametricbem2d::PanelVector;
  Eigen::Vector2d center; center << 2.5,3.6;
  double radius = 33;
  parametricbem2d::ParametrizedCircularArc parametrization(center,radius,.33*M_PI,.99*M_PI);
  unsigned int N = 10;
  PanelVector components = parametrization.split(N);
  for (unsigned int i = 0 ; i < N-1 ; ++i) {
    //std::cout << "i: " << i << std::endl;
    EXPECT_NEAR((components[i]->operator()(1) -
                 components[i+1]->operator()(-1)).norm(),0,eps);
  }
}

TEST(Split,ParametrizedFourierSum) {
  Eigen::MatrixXd a(2,1); //cosine coefficients ; N = 1
  Eigen::MatrixXd b(2,1); //sine coefficients ; N = 1
  a << 1.,
       0.;
  b << 0.,
       1.;
  parametricbem2d::ParametrizedFourierSum parametrization(a,b);
  unsigned int N = 10;
  using PanelVector = parametricbem2d::PanelVector;
  PanelVector components = parametrization.split(N);
  for (unsigned int i = 0 ; i < N-1 ; ++i)
    EXPECT_NEAR((components[i]->operator()(1) -
                 components[i+1]->operator()(-1)).norm(),0,eps);

}

TEST(ParametrizedMeshTest,MemberFunctions) {
  using PanelVector = parametricbem2d::PanelVector;
  Eigen::Vector2d x1; x1 << 0,0; // Point (0,0)
  Eigen::Vector2d x2; x2 << 1,0; // Point (1,0)
  Eigen::Vector2d x3; x3 << 1,1; // Point (1,1)
  Eigen::Vector2d x4; x4 << 0,1; // Point (0,1)
  parametricbem2d::ParametrizedLine line1(x1,x2);
  parametricbem2d::ParametrizedLine line2(x2,x3);
  parametricbem2d::ParametrizedLine line3(x3,x4);
  parametricbem2d::ParametrizedLine line4(x4,x1);
  PanelVector line1panels = line1.split(2);
  PanelVector line2panels = line2.split(2);
  PanelVector line3panels = line3.split(2);
  PanelVector line4panels = line4.split(2);
  PanelVector panels;
  panels.insert(panels.end(),line1panels.begin(),line1panels.end());
  panels.insert(panels.end(),line2panels.begin(),line2panels.end());
  panels.insert(panels.end(),line3panels.begin(),line3panels.end());
  panels.insert(panels.end(),line4panels.begin(),line4panels.end());
  parametricbem2d::ParametrizedMesh mesh(panels);
  Eigen::Vector2d vertex2 = mesh.getVertex(1);
  EXPECT_EQ(mesh.getNumPanels(),8);
  EXPECT_NEAR(vertex2(0),0.5,eps);
  EXPECT_NEAR(vertex2(1),0,eps);
}

TEST(SingleLayer,PanelOrientedAssembly0) {
  using PanelVector = parametricbem2d::PanelVector;
  Eigen::Vector2d x1; x1 << 0,0; // Point (0,0)
  Eigen::Vector2d x2; x2 << 1,0; // Point (1,0)
  Eigen::Vector2d x3; x3 << 1,1; // Point (1,1)
  Eigen::Vector2d x4; x4 << 0,1; // Point (0,1)
  parametricbem2d::ParametrizedLine line1(x1,x2);
  parametricbem2d::ParametrizedLine line2(x2,x3);
  parametricbem2d::ParametrizedLine line3(x3,x4);
  parametricbem2d::ParametrizedLine line4(x4,x1);
  PanelVector line1panels = line1.split(1);
  PanelVector line2panels = line2.split(2);
  PanelVector line3panels = line3.split(3);
  PanelVector line4panels = line4.split(4);
  PanelVector panels;
  panels.insert(panels.end(),line1panels.begin(),line1panels.end());
  panels.insert(panels.end(),line2panels.begin(),line2panels.end());
  panels.insert(panels.end(),line3panels.begin(),line3panels.end());
  panels.insert(panels.end(),line4panels.begin(),line4panels.end());
  parametricbem2d::ParametrizedMesh mesh(panels);
  parametricbem2d::DiscontinuousSpace<0> space;
  Eigen::MatrixXd galerkin = SingleLayerMatrix(mesh,space,32);
  unsigned int numpanels = mesh.getNumPanels();
  EXPECT_EQ(numpanels,galerkin.cols());
  EXPECT_EQ(numpanels,galerkin.rows());
}

TEST(SingleLayer,PanelOrientedAssembly1) {
  using PanelVector = parametricbem2d::PanelVector;
  Eigen::Vector2d x1; x1 << 0,0; // Point (0,0)
  Eigen::Vector2d x2; x2 << 1,0; // Point (1,0)
  Eigen::Vector2d x3; x3 << 1,1; // Point (1,1)
  Eigen::Vector2d x4; x4 << 0,1; // Point (0,1)
  parametricbem2d::ParametrizedLine line1(x1,x2);
  parametricbem2d::ParametrizedLine line2(x2,x3);
  parametricbem2d::ParametrizedLine line3(x3,x4);
  parametricbem2d::ParametrizedLine line4(x4,x1);
  PanelVector line1panels = line1.split(1);
  PanelVector line2panels = line2.split(1);
  PanelVector line3panels = line3.split(1);
  PanelVector line4panels = line4.split(1);
  PanelVector panels;
  panels.insert(panels.end(),line1panels.begin(),line1panels.end());
  panels.insert(panels.end(),line2panels.begin(),line2panels.end());
  panels.insert(panels.end(),line3panels.begin(),line3panels.end());
  panels.insert(panels.end(),line4panels.begin(),line4panels.end());
  parametricbem2d::ParametrizedMesh mesh(panels);
  parametricbem2d::DiscontinuousSpace<1> space;
  Eigen::MatrixXd galerkin = SingleLayerMatrix(mesh,space,32);
  unsigned int numpanels = mesh.getNumPanels();
  EXPECT_EQ(8,galerkin.cols());
  EXPECT_EQ(8,galerkin.rows());
}

TEST(SingleLayer,CppHilbertComparison) {
  using PanelVector = parametricbem2d::PanelVector;
  Eigen::Vector2d x1; x1 << 0,0; // Point (0,0)
  Eigen::Vector2d x2; x2 << 1,0; // Point (1,0)
  Eigen::Vector2d x3; x3 << 1,1; // Point (1,1)
  Eigen::Vector2d x4; x4 << 0,1; // Point (0,1)
  parametricbem2d::ParametrizedLine line1(x1,x2);
  parametricbem2d::ParametrizedLine line2(x2,x3);
  parametricbem2d::ParametrizedLine line3(x3,x4);
  parametricbem2d::ParametrizedLine line4(x4,x1);
  PanelVector line1panels = line1.split(1);
  PanelVector line2panels = line2.split(1);
  PanelVector line3panels = line3.split(1);
  PanelVector line4panels = line4.split(1);
  PanelVector panels;
  panels.insert(panels.end(),line1panels.begin(),line1panels.end());
  panels.insert(panels.end(),line2panels.begin(),line2panels.end());
  panels.insert(panels.end(),line3panels.begin(),line3panels.end());
  panels.insert(panels.end(),line4panels.begin(),line4panels.end());
  parametricbem2d::ParametrizedMesh parametrizedmesh(panels);
  parametricbem2d::DiscontinuousSpace<0> space;
  Eigen::MatrixXd galerkinnew = SingleLayerMatrix(parametrizedmesh,space,32);
  Eigen::MatrixXd coords(4,2); coords << x1,x2,x3,x4;
  Eigen::Matrix<int,4,2> elems;
  elems << 0,1,
           1,2,
           2,3,
           3,0;
  BoundaryMesh boundarymesh(coords,elems);
  Eigen::MatrixXd galerkinold;
  computeV(galerkinold,boundarymesh,0);
  EXPECT_NEAR((galerkinold-galerkinnew).norm(),0,eps);
}

int main(int argc, char **argv) {
  srand(time(NULL));
  // run tests
  ::testing::InitGoogleTest(&argc,argv);
  return RUN_ALL_TESTS();
}
