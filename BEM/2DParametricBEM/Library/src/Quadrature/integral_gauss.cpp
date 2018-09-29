/**
 * \file integral_gauss.cpp
 * \brief This file declares the functions to evaluate the entries of
 *        Galerkin matrices using the transformations given in section
 *        1.4.3.4 in the Lecture Notes for Advanced Numerical Methods
 *        for CSE.
 *
 * This File is a part of the 2D-Parametric BEM package
 */
#include "integral_gauss.hpp"

#include <Eigen/Dense>
#include "genLaguerreRule.hpp"
#include "gaussQuadrature.h"
#include "logweight_quadrature.hpp"

namespace parametricbem2d {
  template <typename T>
  double ComputeIntegral(T integrand,double a,double b,unsigned int N) {
    const double* weights = getGaussWeights(N);
    const double* points = getGaussPoints(N);
    double integral = 0.;
    for ( unsigned int i = 0; i<N ; ++i ) {
      double x = (b-a)/2.*points[i]+(b+a)/2.;
      integral += weights[i]*integrand(x);
    }
    integral *= (b-a)/2;
    return integral;
  }

  template <typename T>
  double ComputeLogIntegral(T integrand,double a,unsigned int N) {
    QuadRule logweightQR = getLogWeightQR(a,N);
    double integral = 0.;
    for ( unsigned int i = 0; i<N ; ++i ) {
      double x = logweightQR.x(i);
      integral += logweightQR.w(i)*integrand(x);
    }
    return integral;
  }
} // namespace parametricbem2d
