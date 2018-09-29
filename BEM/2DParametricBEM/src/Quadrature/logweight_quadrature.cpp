#include "logweight_quadrature.hpp"

#include <Eigen/Dense>
#include "genLaguerreRule.hpp"

static const int KIND = 5;

QuadRule getLogWeightQR(double a, int n){
  // get Gauss-Laguerre quadrature points and weights
  double *w = new double[n];
  double *x = new double[n];
  cgqf(n, KIND, 0, 0, 0, 1, x, w);
  // create new Quadrature rule
  QuadRule logWeightQR;
  logWeightQR.dim = 1;
  logWeightQR.n = n;
  // create matrix and vector for its points and weights
  Eigen::MatrixXd points(n,1);
  Eigen::VectorXd weights(n);
  // Generate them applying the required transformation
  for(int i=0; i<n; i++){
    // quadrature point becomes: $q = a e^{-x}$
    points(i,0) = a*exp(-x[i]);
    // weight becomes $w = a(log(a)-x)w$
    weights(i) = a*(log(a)-x[i])*w[i];
  }
  logWeightQR.x = points;
  logWeightQR.w = weights;
  return logWeightQR;
}
