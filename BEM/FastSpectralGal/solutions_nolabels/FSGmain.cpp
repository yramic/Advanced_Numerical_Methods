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

#include "PeriodicTrapezoidalQR.hpp"



/* @brief Compute matrix A-M using analytic expression.
 * \param[in] N Discretization parameter indicating number of basis functions.
 */
MatrixXd computeM(int N){
  Eigen::MatrixXd M(N,N);
  M.setZero();

// Matrix assembly
for (unsigned int i = 0 ; i < N ; ++i) {
  for (unsigned int j = 0 ; j < N ; ++j) {
    if (i==j)
      M(i,j) = M_PI/4./(i+1);
  }
}
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
  // Get periodic Trapezoidal Quadrature points and weight
  VectorXd TR_points(2*N);  double TR_w;
  std::tie(TR_points,TR_w) = PeriodicTrapRule(2*N);
  // Fill vector entries
  for (unsigned int i = 0 ; i < N ; ++i) {
    // Evaluate integral by using the quadrature
    for (unsigned int qp = 0 ; qp < 2*N ; ++qp) {
      double a = TR_points(qp);
      RHS(i) += 0.5 * TR_w * g(cos(a/2))*cos((i+1)*a/2.);
    }// end iteration over quadrature points
  }

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

  std::cout << " Assemble M " << std::endl;
  Eigen::MatrixXd M = computeM(N);
  // Build RHS
  std::cout << " Assemble RHS" << std::endl;
  VectorXd RHS = computeG(g, N);
  // Use direct solver
  std::cout << " Solve " << std::endl;
  MatrixXd LHS = M.eval();
  VectorXd sol = LHS.lu().solve(RHS);
  std::cout << " Done " << std::endl;

  return sol;
}


//----------------------------------------------------------------------------
/* @brief Build and solve boundary integral equation V rho = g ignoring first
 *        row of the system (A has empty row for on the Disk).
 * \param[in] gamma Function that takes a double and returns a 2d vector
 *                  corresponding to the parametrized curve.
 * \param[in] g Right hand side function (takes 2d points and returns a double).
 * \param[in] N Discretization parameter indicating number of basis functions.
 */
template <typename PARAM, typename FUNC>
VectorXd solveBIEonDisk(const PARAM& gamma, const FUNC& g, int N){
  // In order to compute A we do (A-M)+M
  std::cout << " Assemble A " << std::endl;
  MatrixXd AmM = computeAminusM(N);
  MatrixXd M   = computeM(gamma, N);
  // Build RHS
  std::cout << " Assemble RHS" << std::endl;
  VectorXd RHS = computeG(gamma, g, N);
  std::cout << " Solve " << std::endl;
  MatrixXd LHS = (AmM + M).block(1,1,2*N,2*N);
  VectorXd sol = LHS.lu().solve(RHS.segment(1,2*N));
  std::cout << " Done " << std::endl;
  return sol;
}


//----------------------------------------------------------------------------
/* @brief Reconstruct function UN from its coefficients and evaluate it at t.
 * \param[in] coeffs coefficients of UN
 * \param[in] t point in [0,2Pi]
 * \param[in] gammaprime Function that takes a double and returns a 2d vector
 *                       corresponding to the derivative of the curve's
 *                       parametrization.
 */
template <typename PARAMDER>
double reconstructRho(const VectorXd& coeffs, double t){
  int N = coeffs.rows();
  double rho_N = 0.;
  for (unsigned int i = 0 ; i < N ; ++i) {
    rho_N += coeffs(i) * cos((i+1) * acos(t));
  }
  return rho_N/std::sqrt(1-t*t);
}


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
  int N = coeffs.rows();
  // Get quadrature points and weight
  VectorXd TR_points(2*N);  double TR_w;
  std::tie(TR_points,TR_w) = PeriodicTrapRule(2*N);
  for(int qp=0; qp<2*N; qp++){
    auto z = TR_points(qp);
    // evaluate
    double rho = reconstructRho(coeffs, z);
    norm += TR_w*rho*rho;
  }

  return std::sqrt(res);
}


int main() {
  unsigned int N = 10;
  Eigen::MatrixXd M = computeM(N);
  std::cout << "M is \n" << M << std::endl;
  return 0;

}
