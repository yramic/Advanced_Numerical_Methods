/**
 * \file integral_gauss.hpp
 * \brief This file declares the functions to evaluate the entries of
 *        Galerkin matrices using the transformations given in section
 *        1.4.3.4 in the Lecture Notes for Advanced Numerical Methods
 *        for CSE.
 *
 * This File is a part of the 2D-Parametric BEM package
 */

#ifndef INTEGRALGAUSSHPP
#define INTEGRALGAUSSHPP

#include <Eigen/Dense>

namespace parametricbem2d {
  template <typename T>
  double ComputeIntegral(T integrand,double a,double b,unsigned int N);

  template <typename T>
  double ComputeLogIntegral(T integrand,double a,unsigned int N);
} // namespace parametricbem2d

#endif //INTEGRALGAUSSHPP
