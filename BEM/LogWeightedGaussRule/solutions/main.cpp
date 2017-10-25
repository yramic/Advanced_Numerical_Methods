#include <iostream>
#include <fstream>
#include <istream> 
#include <cmath>
#include <Eigen/Dense>

#include "genLaguerreRule.hpp"

// Quadrature rule struct
struct QuadRule {
  std::size_t     dim; // dimension of space
  std::size_t     n;   // number of nodes/weights
  Eigen::MatrixXd x;   // quadrature nodes (columns of a matrix with dim rows)
  Eigen::VectorXd w;   // vector of quadrature weights
};

static const int KIND = 5;

//------------------------------------------------------------------------------
/* SAM_LISTING_BEGIN_0 */
bool testGaussLaguerre(int N){
  // get Gauss-Laguerre quadrature points and weights
  double *w = new double[N];
  double *x = new double[N];
  cgqf(N, KIND, 0, 0, 0, 1, x, w);

  // initialize exact value of integral
  double exval = 1;
  // Start test for different k=1..N
  for(int k=0; k<2*N; k++){
    // compute exact value for $\int_0^{\infty} e^{-t} t^k dt$
    if(k>0)
      exval *= k;
    // compute integral using quadrature
    double I = 0.;
    for(int i=0; i<N; i++){
      I += std::pow(x[i], k)*w[i];
    }
    if(fabs(I-exval)>1e-10*fabs(exval)){ // if not exact, stop
      return false;
    }
  }
  return true;
}
/* SAM_LISTING_END_0 */


//------------------------------------------------------------------------------
/* SAM_LISTING_BEGIN_1 */
Eigen::VectorXd testGaussLaguerreConvergence(int N = 20){
  Eigen::VectorXd error(N);
  for(int k=1; k<=N; k++){
    // get Gauss-Laguerre quadrature points and weights
    double *w = new double[k];
    double *x = new double[k];
    cgqf(k, KIND, 0, 0, 0, 1, x, w);
    // compute integral using quadrature
    double I = 0.;
    for(int i=0; i<k; i++){
      I += sin(x[i])*w[i];
    }
    // compute error
    error(k-1) = fabs(I - 0.5);
  }
  return error;
}
/* SAM_LISTING_END_1 */

//------------------------------------------------------------------------------
/* SAM_LISTING_BEGIN_2 */
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
/* SAM_LISTING_END_2 */


//------------------------------------------------------------------------------
/* SAM_LISTING_BEGIN_3 */
Eigen::VectorXd testLogWeightQRConvergence(int n){
  Eigen::VectorXd error(n);
  for(int k=1; k<=n; k++){
    // get quadrature points and weights
    QuadRule LWQR = getLogWeightQR(1, k);
    // compute integral using quadrature
    double I = 0.;
    for(int i=0; i<k; i++){
      I += sin(LWQR.x(i))*LWQR.w(i);
    }
    // compute error
    error(k-1) = fabs(I + 0.239811742000564);
  }
  return error;
}
/* SAM_LISTING_END_3 */


//------------------------------------------------------------------------------
/* SAM_LISTING_BEGIN_4 */
int testLogWeightQR(int N){
  // initialize exact value
  double exval;
  // get quadrature points and weights
  QuadRule LWQR = getLogWeightQR(1, N);
  for(int k=1; k<2*N; k++){
    // compute exact value for $\int_0^1 log(t) t^k dt$
    exval = -1./std::pow(k+1,2);
    // compute integral using quadrature
    double I = 0.;
    for(int i=0; i<N; i++){
      I += std::pow(LWQR.x(i), k)*LWQR.w(i);
    }
    if(fabs(I-exval)>1e-10*fabs(exval)){ // if not exact, stop
      return k;
    }
  }
  return N;
}
/* SAM_LISTING_END_4 */


//------------------------------------------------------------------------------
int main() {

  //-------------------------------------------------------------
  /* SAM_LISTING_BEGIN_5 */
  std::cout << "TESTING GAUSS-LAGUERRE QUADRATURE" << std::endl;
  for(int j=1; j<30; j++){
    bool pass = testGaussLaguerre(j);
    if(pass){
      std::cout << "Test passed for N = " << j << std::endl;
    }
    else{
      std::cout << "Test for N = " << j << " failed " << std::endl;
    }
  }
  /* SAM_LISTING_END_5 */
  std::cout << std::endl;

  //-------------------------------------------------------------
  /* SAM_LISTING_BEGIN_6 */
  std::cout << "TESTING QUADRATURES' CONVERGENCE" << std::endl
	    << std::endl;
  Eigen::VectorXd errorGL = testGaussLaguerreConvergence(30);
  Eigen::VectorXd errorLW = testLogWeightQRConvergence(30);
  Eigen::VectorXi N = Eigen::VectorXi::LinSpaced(30, 1, 30);
  //output for plot
  std::ofstream out_errorGL("GLQR_errors.txt");
  out_errorGL << std::setprecision(18) << errorGL; 
  out_errorGL.close( );
  std::ofstream out_errorLW("WLQR_errors.txt");
  out_errorLW << std::setprecision(18) << errorLW; 
  out_errorLW.close( );
  std::ofstream out_N("QR_N.txt");
  out_N << N; 
  out_N.close( );
  /* SAM_LISTING_END_6 */

  //-------------------------------------------------------------
  /* SAM_LISTING_BEGIN_7 */
  std::cout << "TESTING LOG-WEIGHT QUADRATURE" << std::endl;
  for(int j=1; j<30; j++){
    int orderfail = testLogWeightQR(j);
    if(orderfail >= j){
      std::cout << "Test passed for N = " << j  << std::endl;
    }
    else{
      std::cout << "Test for N = " << j << " failed at " << orderfail
		<< std::endl;
    }
  }
  /* SAM_LISTING_END_7 */
  std::cout << std::endl;
  
  return 0;

}
