/**
 * @file mfmg.cpp
 * @brief ADVNCSE homework MFMG code
 * @author Bob Schreiner, Jörg Nick
 * @date October 2023
 * @copyright Developed at SAM, ETH Zurich
 */

#include "mfmg.h"

#include <Eigen/Sparse>
#include <unsupported/Eigen/KroneckerProduct>
namespace MFMG {

/* SAM_LISTING_BEGIN_1 */
GridFunction DirichletBVPMultiGridSolver::applyGridOperator(
    unsigned int level, const GridFunction &mu) const {
  // Assume the gridfunction $mu$ vanishes at the boundary
  const unsigned int M = std::pow(2, level) - 1;
  const double h = 1.0 / (M + 1);
  GridFunction eta = Eigen::MatrixXd::Zero(M + 2, M + 2);
  // ************************************************************
  // PROBLEM 4-2C:
  for (unsigned int i{1}; i < M+1; ++i) {
    for (unsigned int j{1}; j < M+1; ++j) {
      eta(i, j) = (4.0 + c_ * std::pow(h,2.0)) * mu(i,j) - 
                  mu(i-1,j) - mu(i,j-1) - mu(i+1,j) -mu(i,j+1);
    }
  } 
  // ************************************************************
  return eta;
}
/* SAM_LISTING_END_1 */

/* SAM_LISTING_BEGIN_2 */
GridFunction DirichletBVPMultiGridSolver::directSolve(
    unsigned int level, const GridFunction &phi) const {
  const unsigned int M = std::pow(2, level) - 1;
  const double h = 1.0 / (M + 1);
  const Eigen::MatrixXd phi_inner = phi.block(1, 1, M, M);
  const Eigen::VectorXd phi_v =
      Eigen::MatrixXd::Map(phi_inner.data(), M * M, 1);
  Eigen::VectorXd sol(M * M);
  Eigen::SparseMatrix<double> A(M * M, M * M);
  Eigen::SparseMatrix<double> A1(M * M, M * M);
  Eigen::SparseMatrix<double> A2(M * M, M * M);
  Eigen::SparseMatrix<double> T(M, M);
  A.reserve(Eigen::VectorXi::Constant(A.cols(), 5));
  A1.reserve(Eigen::VectorXi::Constant(A1.cols(), 5));
  A2.reserve(Eigen::VectorXi::Constant(A2.cols(), 5));
  T.reserve(Eigen::VectorXi::Constant(T.cols(), 3));
  // ************************************************************
  // PROBLEM 4-2D:
  // First we need to define T properly! This cone be done with
  // T.insert()
  // Define the first row:
  T.insert(0,0) = 2.0 + 0.5*c_*std::pow(h,2.0);
  T.insert(0,1) = -1;
  // Define the last row:
  T.insert(M-1,M-1) = 2.0 + 0.5*c_*std::pow(h,2.0);
  T.insert(M-1,M-2) = -1;
  // Now we can run a Loop and define the tridiagonal elements of T:
  for (unsigned int i{1}; i < M-1; ++i) {
    T.insert(i,i) = 2.0 + 0.5*c_*std::pow(h,2.0);
    T.insert(i,i-1) = -1;
    T.insert(i,i+1) = -1;
  }

  // Next we want to define the Matrix A, which can be done with
  // the Kronecker Product!

  Eigen::KroneckerProductSparse kron_1(T, 
        Eigen::MatrixXd::Identity(T.rows(), T.cols()));
  Eigen::KroneckerProductSparse kron_2(
        Eigen::MatrixXd::Identity(T.rows(), T.cols()), T);
  
  // Now we can Assign these results to A1 and A2:
  kron_1.evalTo(A1);
  kron_2.evalTo(A2);

  A = A1 + A2;

  // Now we have everything we need for the Sparse Solver and we
  // can finally compute the result directly on the coarsest grid:
  Eigen::SparseLU<Eigen::SparseMatrix<double>> sparse_solver;
  sparse_solver.analyzePattern(A);
  sparse_solver.factorize(A);
  assert(sparse_solver.info() == Eigen::Success);
  sol = sparse_solver.solve(phi_v);
  // ************************************************************
  GridFunction sol_g = Eigen::MatrixXd::Zero(M + 2, M + 2);
  sol_g.block(1, 1, M, M) = Eigen::MatrixXd::Map(sol.data(), M, M);
  return sol_g;
}
/* SAM_LISTING_END_2 */

/* SAM_LISTING_BEGIN_3 */
void DirichletBVPMultiGridSolver::sweepGaussSeidel(
    unsigned int level, GridFunction &mu, const GridFunction &phi) const {
  const unsigned int M = std::pow(2, level) - 1;
  const double h = 1.0 / (M + 1);
  // ************************************************************
  // PROBLEM 4-2e:
  for (unsigned int i{1}; i < M+1; ++i) {
    for (unsigned int j{1}; j < M+1; ++j) {
      mu(i,j) = (phi(i,j) + (mu(i-1,j) + mu(i,j-1) + mu(i+1,j) + 
                mu(i,j+1))) / (4 + c_ * std::pow(h, 2.0));
    }
  }
  // ************************************************************
}
/* SAM_LISTING_END_3 */

/* SAM_LISTING_BEGIN_4 */
GridFunction DirichletBVPMultiGridSolver::residual(
    unsigned int level, const GridFunction &mu, const GridFunction &phi) const {
  const unsigned int M = std::pow(2, level) - 1;
  GridFunction res(M + 2, M + 2);
  // ************************************************************
  // PROBLEM 4-2F:
  res = (phi - applyGridOperator(level, mu));
  // ************************************************************
  return res;
}

/* SAM_LISTING_END_4 */

/* SAM_LISTING_BEGIN_5 */
GridFunction DirichletBVPMultiGridSolver::prolongate(
    unsigned int level, const GridFunction &gamma) const {
  const unsigned int M = std::pow(2, level) - 1;
  GridFunction zeta((M + 2), (M + 2));
  zeta.setZero();
  // ************************************************************
  // PROBLEM 4-2G:
  // A nested for loop is requried:
  for (unsigned int i{1}; i < M+1; ++i) {
    for (unsigned int j{1}; j < M+1; ++j) {
      // Now I want to check the 4 cases:
      // Case 1: if i,j == even:
      if ((i % 2 == 0) and (j % 2 == 0)) 
        zeta(i,j) = gamma(i/2, j/2);
      // Case 2: if i == even, j == odd:
      if ((i % 2 == 0) and (j % 2 != 0))
        zeta(i,j) = 0.5 * (gamma(i/2, (j+1)/2) + gamma(i/2, (j-1)/2));
      // Case 3: if i == odd, j == even
      if((i % 2 != 0) and (j % 2 == 0))
        zeta(i,j) = 0.5 * (gamma((i-1)/2, j/2) + gamma((i+1)/2, j/2));
      // Case 4: if i,j == odd:
      if ((i % 2 != 0) and (j % 2 != 0)) {
        zeta(i,j) = 0.25 * (gamma((i-1)/2, (j-1)/2) + gamma((i-1)/2, (j+1)/2) +
                            gamma((i+1)/2, (j-1)/2) + gamma((i+1)/2, (j+1)/2));
      }
    }
  }
  // ************************************************************
  return zeta;
}

/* SAM_LISTING_END_5 */

/* SAM_LISTING_BEGIN_6 */
GridFunction DirichletBVPMultiGridSolver::restrict(
    unsigned int level, const GridFunction &rho) const {
  const unsigned int M = std::pow(2, level - 1) - 1;
  GridFunction sigma((M + 2), (M + 2));
  sigma.setZero();
  // ************************************************************
  // PROBLEM 4-2H:
  for (unsigned int i{1}; i < M+1; ++i) {
    for (unsigned int j{1}; j < M+1; ++j) {
      sigma(i,j) = rho(2*i, 2*j);
      sigma(i,j) += 0.5 * (rho(2*i,2*j-1) + rho(2*i, 2*j+1) + 
                           rho(2*i-1, 2*j) + rho(2*i+1, 2*j));
      sigma(i,j) += 0.25 * (rho(2*i-1, 2*j-1) + rho(2*i-1, 2*j+1) + 
                            rho(2*i+1, 2*j-1) + rho(2*i+1, 2*j+1));
    }
  }
  // ************************************************************
  return sigma;
}
/* SAM_LISTING_END_6 */

/* SAM_LISTING_BEGIN_7 */
GridFunction DirichletBVPMultiGridSolver::multigridIteration(
    const GridFunction &mu, const GridFunction &phi, unsigned int L0) const {
  std::vector<GridFunction> mu_vec;
  std::vector<GridFunction> phi_vec;
  std::vector<GridFunction> residual_vec;
  mu_vec.push_back(mu);
  const unsigned int M = std::pow(2, L0) - 1;
  const double h = 1.0 / (M + 1);
  phi_vec.push_back(phi);

  // ************************************************************
  // PROBLEM 4-2I:
  unsigned counter {0};
  // We go from level l to level 0 where we use a direct solver
  // Hence, the forloop will go actually backward starting from
  // l and ending at level 1, because level 0 needs to be solved
  // directly as mentioned:
  for (unsigned int l {L_}; l > L0; --l) {
    // Pre-Smoothing:
    sweepGaussSeidel(l, mu_vec.back(), phi_vec.back());
    // Compute the residual:
    residual_vec.push_back(
      residual(l, mu_vec.back(), phi_vec.back())
    );
    // Restirct resiudal:
    phi_vec.push_back(
      restrict(l, residual_vec.back())
    );
    // Add zero grid function to mu:
    GridFunction zero {phi_vec.back()};
    zero.setZero();
    mu_vec.push_back(zero);
    ++counter;
  }
  // At Level L0 we do a direct solver:
  mu_vec.push_back(directSolve(L0, phi_vec.back()));
  // Now Prolongate and do a post smoothening:
  for (unsigned int l{L0}; l < L_; ++l) {
    // Prolongate mu:
    mu_vec[counter - 1] += prolongate(l+1, mu_vec[counter]);
    // Post Smoothening:
    sweepGaussSeidel(l+1, mu_vec[counter-1], phi_vec[counter-1]);
    --counter;
  }
  // ************************************************************
  return mu_vec[0];
}
/* SAM_LISTING_END_7 */

/* SAM_LISTING_BEGIN_8 */
double estimateMGConvergenceRate(double c, unsigned int L, double tol) {
  unsigned M = std::pow(2, L) - 1;
  double h = 1. / (M + 1);
  double lambda_new = 1.;
  double lambda_old = 0.;
  DirichletBVPMultiGridSolver solver(c, L);
  GridFunction v = Eigen::MatrixXd::Constant(M + 2, M + 2, 1);
  v.col(0) = Eigen::VectorXd::Zero(M + 2);
  v.col(M + 1) = Eigen::VectorXd::Zero(M + 2);
  v.row(0) = Eigen::VectorXd::Zero(M + 2);
  v.row(M + 1) = Eigen::VectorXd::Zero(M + 2);

  GridFunction v_old;
  GridFunction Av;
  // ************************************************************
  // Code to be supplemented
  // ************************************************************
  return lambda_new;
}
/* SAM_LISTING_END_8 */

/* SAM_LISTING_BEGIN_9 */
void tabulateMGConvergenceRate() {
  std::vector<double> c_vec({-40, -20, -10, -1, 0, 1, 10, 20, 40});
  std::vector<unsigned> L_vec({6, 7, 8, 9, 10, 11});
  // ************************************************************
  // Code to be supplemented
  // ************************************************************
}
/* SAM_LISTING_END_9 */

}  // namespace MFMG
