#include <iostream>
#include <cmath>
#include <Eigen/Dense>


/* @brief Compute Periodic Trapezoidal Rule over unit circle
 * \param[in] n Number of quadrature points
 */
std::pair<Eigen::MatrixXd,double>
PeriodicTrapRule(int N) {

  double w1 = 2*M_PI/N;
 
  Eigen::MatrixXd xq(N,2); // quadrature points

  for(int i=0; i<N; i++){
    xq.row(i) << cos(2*M_PI*j/N), sin(2*M_PI*j/N);
  }

  // return
  return std::make_pair(xq,wq);

}
