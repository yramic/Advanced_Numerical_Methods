////
//// Copyright (C) 2017 SAM (D-MATH) @ ETH Zurich
//// Author(s): curzuato < >
//// Contributors:  dcasati
//// This file is part of the AdvNumCSE repository.
////
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <Eigen/SparseLU>
#include <cassert>
#include <cmath>
#include <ctime>
#include <iomanip>
#include <iostream>
#include <unsupported/Eigen/FFT>

using namespace Eigen;
using namespace std;

/* @brief Galerkin matrix for FE discretisation of $-\nabla$
 * for a uniform triangulation on an equilateral triangle
 * \param l refinement level
 * \\return sparse Galerkin matrix
 */
SparseMatrix<double> genGalerkinMat(const int l) {
  // TODO: Implement the Galerkin matrix generator
}

/* @brief Prolongation matrix $\VP_{\ell-1, \ell}$
 * \param l refinement level
 * \\return sparse Prolongation matrix
 */
SparseMatrix<double> genProlongationMat(const int l) {
  // TODO: Implement the prolongation matrix generator
}

/* @brief Multi-grid method
 * to solve sparse linear system of equations
 * \param phi right-hand-side vector, $\IR^{N}$
 * \param mu initial guess and output solution, $\IR^{N \times N}$
 * \param max_n_steps maximum number of steps
 * \param TOL error tolerance for termination criteria
 */
void multiGridIter(const VectorXd &phi, VectorXd &mu, int l, int max_n_steps,
                   double TOL = 1E-04) {
  // TODO: Implement multi-grid iteration
}

/* @brief Computes the residual for the Multi-grid method with max_n_steps=1
 * \param L mesh level
 * \param max_itr maximum number of multi-grid iterations
 * \param TOL error tolerance for termination criteria
 */
VectorXd test_multigrid_residual(const int L, int max_itr,
                                 double TOL = 1.0E-08) {
  VectorXd rho_norm(max_itr);

  unsigned n = std::pow(2, L);
  unsigned N = 0.5 * (n + 2) * (n + 1);  // Number of nodes

  SparseMatrix<double> A = genGalerkinMat(L);
  VectorXd phi = VectorXd::Zero(N);
  VectorXd mu = VectorXd::Random(N);

  unsigned int count = 0;
  for (int i = 0; i < max_itr; i++) {
    multiGridIter(phi, mu, L, 1);
    rho_norm(i) = (phi - A * mu).norm();

    count++;
    if (rho_norm(i) <= TOL) break;
  }

  return rho_norm.head(count);
}

/* @brief Computes the asymptotic convergence rates of the Multi-grid method
 * by measuring the quotient of the residuals from successive iterations.
 */
void test_multigrid_convergence() {
  int max_itr = 100;
  cout << "L\trate" << endl;
  for (int L = 2; L <= 8; L++) {
    VectorXd rho_norm = test_multigrid_residual(L, max_itr);
    VectorXd rates(rho_norm.size() - 1);
    for (int i = 0; i < rates.size(); i++)
      rates(i) = rho_norm(i + 1) / rho_norm(i);
    cout << L << "\t" << rates.tail(1) << endl;
  }
}

int main() {
  int l = 2;
  double TOL = 1E-04;

  SparseMatrix<double> A = genGalerkinMat(l);
  cout << "\n\nGalerkin matrix for l=" << l << ":\n\n" << A << endl;

  SparseMatrix<double> P = genProlongationMat(l);
  cout << "\n\nProlongation matrix for l=" << l << ":\n\n" << P << endl;

  unsigned n = std::pow(2, l);
  unsigned N = 0.5 * (n + 2) * (n + 1);  // Number of nodes
  VectorXd phi = VectorXd::Zero(N);
  VectorXd mu = VectorXd::Random(N);

  int max_nsteps = 1000;
  multiGridIter(phi, mu, l, max_nsteps, TOL);
  //cout << "\n" << mu.transpose() << endl;
}

// End of file
