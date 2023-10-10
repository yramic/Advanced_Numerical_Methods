////
//// Copyright (C) 2017 SAM (D-MATH) @ ETH Zurich
//// Author(s): curzuato < >
//// Contributors:  dcasati
//// This file is part of the AdvNumCSE repository.
////
#include <Eigen/Dense>
#include <cassert>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <istream>

using namespace Eigen;

/* @brief Compute matrix M using analytic expression.
 *
 * \param[in] N Discretization parameter indicating number of basis functions.
 */
MatrixXd computeM(unsigned int N) {
  Eigen::MatrixXd M(N, N);
  M.setZero();

  // TODO: Compute M
  return M;
}

/* @brief Compute right hand side using Chebyshev Gauss Quadrature
 *
 * \tparam FUNC Type for the Dirichlet B.C. supporting : g(const double&)
 * \param[in] g Dirichlet boundary condition; Signature : double(const double&)
 * \param[in] N Discretization parameter indicating number of basis functions.
 */
template <typename FUNC>
VectorXd computeG(const FUNC &g, int N) {
  // Initialize right hand side vector
  VectorXd RHS(N);
  RHS.setZero();
  // TODO: Compute RHS
  return RHS;
}

/* @brief Build and solve boundary integral equation V rho = g
 *
 * \tparam FUNC Type for the Dirichlet B.C. supporting : g(const double&)
 * \param[in] g Dirichlet boundary condition; Signature : double(const double&)
 * \param[in] N Discretization parameter indicating number of basis functions.
 */
template <typename FUNC>
VectorXd solveBIE(const FUNC &g, int N) {
  // TODO: Build BIE system and solve it
}

/* @brief Compute L2 norm of UN from its coefficients using Gauss Legendre
 *        Quadrature rule
 *
 * \param[in] coeffs coefficients of UN
 */
double L2norm(const VectorXd &coeffs) {
  double norm = 0.;
  // TODO: implement the computation of the L2norm
  return std::sqrt(norm);
}

int main() {
  std::cout << "Test for source term g1(x) = sin(2*Pi*x)" << std::endl;
  std::cout << "N" << std::setw(15) << "L2error" << std::endl;
  std::cout << "############################" << std::endl;

  // std::function object for source term g1 using a lambda expression
  std::function<double(const double &)> g1 = [](const double &x) {
    return sin(2 * M_PI * x);
  };
  // TODO: implement the convergence test
  return 0;
}
