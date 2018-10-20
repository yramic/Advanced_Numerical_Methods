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

/* @brief Reconstruct function UN from its coefficients and evaluate it at t.
 *
 * \param[in] coeffs coefficients of UN
 * \param[in] t evaluation point [-1,1]
 */
double reconstructRho(const VectorXd &coeffs, double t) {
  assert(t >= -1 && t <= 1); // Asserting evaluation is within the domain
  int N = coeffs.rows();
  double rho_N = 0.;
  for (unsigned int i = 0; i < N; ++i)
    // Coefficients start from $T_1(x)$
    rho_N += coeffs(i) * boost::math::chebyshev_t(i + 1, t);
  return rho_N / std::sqrt(1 - t * t);
}

/* @brief Compute L2 norm of UN from its coefficients using Gauss Legendre
 *        Quadrature rule
 *
 * \param[in] coeffs coefficients of UN
 */
double L2norm(const VectorXd &coeffs) {
  double norm = 0.;
  int N = coeffs.rows();
  // Get quadrature points and weight for Gauss Legendre Quadrature
  unsigned int order = 2 * N; // Quadrature order
  Eigen::RowVectorXd weights, points;
  std::tie(points, weights) = gauleg(-1, 1, order);
  // Iterating over quadrature points
  for (int qp = 0; qp < order; qp++) {
    auto z = points(qp);
    // evaluating the function
    double rho = reconstructRho(coeffs, z);
    norm += weights(qp) * rho * rho;
  }
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
    std::cout << N << std::setw(15) << l2error << std::endl;

  }
  return 0;
}
