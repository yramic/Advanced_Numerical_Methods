#include <iostream>
#include <cmath>
#include <Eigen/Dense>


/* @brief Compute Periodic Trapezoidal Rule over I=[0,2 PI]
 * \param[in] n Number of quadrature points
 */
std::pair<Eigen::VectorXd,double>
PeriodicTrapRule(int N) {

  double wq = 2*M_PI/N;
 
  Eigen::VectorXd aux;
  aux.setLinSpaced(N+1,0,1);
  //std::cout << aux.segment(0,N) << std::endl;
  Eigen::VectorXd xq = aux.segment(0,N)*2*M_PI;

  // return
  return std::make_pair(xq,wq);

}
