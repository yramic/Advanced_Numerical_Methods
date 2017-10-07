#include <iostream>
#include <cmath>
#include <Eigen/Dense>


/* @brief Compute Periodic Trapezoidal Rule over unit circle
 * \param[in] n Number of quadrature points
 */
/*
std::pair<Eigen::MatrixXd,double>
PeriodicTrapRule(int N) {

  double wq = 1./N;
 
  Eigen::MatrixXd xq(N,2); // quadrature points

  for(int i=0; i<N; i++){
    xq.row(i) << cos(2*M_PI*i/N), sin(2*M_PI*i/N);
  }

  // return
  return std::make_pair(xq,wq);

}
*/

/* @brief Compute Periodic Trapezoidal Rule over I=[0,1]
 * \param[in] n Number of quadrature points
 */
std::pair<Eigen::VectorXd,double>
PeriodicTrapRule(int N) {

  double wq = 1./N;
 
  Eigen::VectorXd aux;
  aux.setLinSpaced(1,0,N);
  Eigen::VectorXd xq = aux.inverse();

  // return
  return std::make_pair(xq,wq);

}
