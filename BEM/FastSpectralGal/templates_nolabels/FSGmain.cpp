//// 
//// Copyright (C) 2017 SAM (D-MATH) @ ETH Zurich
//// Author(s): curzuato < > 
//// Contributors:  dcasati 
//// This file is part of the AdvNumCSE repository.
////
#include <iostream>
#include <fstream>
#include <istream>
#include <cmath>
#include <Eigen/Dense>

using namespace Eigen;




/* @brief Compute matrix A-M using analytic expression.
 * \param[in] N Discretization parameter indicating number of basis functions.
 */
MatrixXd computeM(int N){
  Eigen::MatrixXd M(N,N);
  M.setZero();

  // TODO: Compute M
  return M;
}

//----------------------------------------------------------------------------
/* @brief Compute right hand side using periodic trapezoidal rule (2N points).
 * \param[in] gamma Function that takes a double and returns a 2d vector
 *                  corresponding to the parametrized curve.
 * \param[in] g Right hand side function (takes 2d points and returns a double).
 * \param[in] N Discretization parameter indicating number of basis functions.
 */
template <typename PARAM, typename FUNC>
VectorXd computeG(const FUNC& g, int N){
  // Initialize right hand side vector
  VectorXd RHS(N);  RHS.setZero();
    // TODO: Compute RHS
  return RHS;
}


//----------------------------------------------------------------------------
/* @brief Build and solve boundary integral equation V rho = g
 * \param[in] gamma Function that takes a double and returns a 2d vector
 *                  corresponding to the parametrized curve.
 * \param[in] g Right hand side function (takes 2d points and returns a double).
 * \param[in] N Discretization parameter indicating number of basis functions.
 */
template <typename PARAM, typename FUNC>
VectorXd solveBIE(const FUNC& g, int N){
    // TODO: Build BIE system and solve it

  return sol;
}


//----------------------------------------------------------------------------


//----------------------------------------------------------------------------
/* @brief Compute L2 norm of UN from its coefficients using periodic trapezoidal
 *        rule (2N points).
 * \param[in] gammaprime Function that takes a double and returns a 2d vector
 *                       corresponding to the derivative of the curve's
 *                       parametrization.
 * \param[in] coeffs coefficients of UN
 */
template <typename PARAMDER>
double L2norm(const VectorXd& coeffs){
  double norm = 0.;
    // TODO: reconstruct the function and compute its L2norm

  return std::sqrt(res);
}


int main() {
  unsigned int N = 10;
  Eigen::MatrixXd M = computeM(N);
  std::cout << "M is \n" << M << std::endl;
  return 0;

}
