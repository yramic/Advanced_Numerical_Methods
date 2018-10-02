/**
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
