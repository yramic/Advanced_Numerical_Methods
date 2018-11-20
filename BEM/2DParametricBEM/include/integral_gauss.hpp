/**
 * This File is a part of the 2D-Parametric BEM package
 */

#ifndef INTEGRALGAUSSHPP
#define INTEGRALGAUSSHPP

#include "gauleg.hpp"
#include "logweight_quadrature.hpp"
#include <limits>

namespace parametricbem2d {
  template <typename T>
  double ComputeIntegral(T integrand,double a,double b,unsigned int N) {
    // Getting quadrature weights and points
    Eigen::RowVectorXd weights,points;
    std::tie(points,weights) = gauleg(a,b,N,std::numeric_limits<double>::epsilon());
    double integral = 0.;
    for ( unsigned int i = 0; i<N ; ++i )
      integral += weights(i)*integrand(points(i));
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

#endif //INTEGRALGAUSSHPP
