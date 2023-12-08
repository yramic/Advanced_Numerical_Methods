/**
 * @file mfmgLshape.cpp
 * @brief ADVNCSE MFMG code
 * @author Bob Schreiner, Jörg Nick
 * @date October 2023
 * @copyright Developed at SAM, ETH Zurich
 */

#include <Eigen/Sparse>
#include <unsupported/Eigen/KroneckerProduct>

#include "mfmg.h"
namespace MFMGLshape {
/* SAM_LISTING_BEGIN_1 */
GridFunction DirichletBVPMultiGridSolver::applyGridOperator(
    unsigned int level, const GridFunction &mu) const {
  // Assume the gridfunction $mu$ vanishes at the boundary
  const unsigned int M = std::pow(2, level) - 1;
  const double h = 1.0 / (M + 1);
  GridFunction eta = Eigen::MatrixXd::Zero(M + 2, M + 2);
  // ************************************************************
  // Code to be supplemented
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
  // Unit matrix eye
  Eigen::SparseMatrix<double> Eye(M, M);
  A.reserve(Eigen::VectorXi::Constant(A.cols(), 5));
  A1.reserve(Eigen::VectorXi::Constant(A1.cols(), 5));
  A2.reserve(Eigen::VectorXi::Constant(A2.cols(), 5));
  T.reserve(Eigen::VectorXi::Constant(T.cols(), 3));
  Eye.reserve(Eigen::VectorXi::Constant(T.cols(), 1));
  T.insert(0, 0) = 2 + (h * h) * 0.5 * c_;
  T.insert(0, 1) = -1;
  T.insert(M - 1, M - 1) = 2. + (h * h) * 0.5 * c_;
  T.insert(M - 1, M - 2) = -1.;
  Eye.insert(0, 0) = 1.;
  Eye.insert(M - 1, M - 1) = 1.;
  for (unsigned int i = 1; i < M - 1; ++i) {
    T.insert(i, i) = 2. + (h * h) * 0.5 * c_;
    T.insert(i, i + 1) = -1.;
    T.insert(i, i - 1) = -1.;
    Eye.insert(i, i) = 1.;
  }
  Eigen::KroneckerProductSparse kron1(Eye, T);
  Eigen::KroneckerProductSparse kron2(T, Eye);
  kron1.evalTo(A1);
  kron2.evalTo(A2);
  A = A1 + A2;
  int counter = 0;
  // middle of the square, from north to south and east to west
  int midpoint = (M + 1) / 2;

  // ************************************************************
  // Code to be supplemented
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
  const double midpoint = (M + 1) / 2;
  // ************************************************************
  // Code to be supplemented
  // ************************************************************
}
/* SAM_LISTING_END_3 */

/* SAM_LISTING_BEGIN_4 */
GridFunction DirichletBVPMultiGridSolver::residual(
    unsigned int level, const GridFunction &mu, const GridFunction &phi) const {
  const unsigned int M = std::pow(2, level) - 1;
  GridFunction res(M + 2, M + 2);
  // ************************************************************
  // Code to be supplemented
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
  // Code to be supplemented
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
  // Code to be supplemented
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
  // Code to be supplemented
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
  int midpoint = (M + 1) / 2;
  GridFunction v = Eigen::MatrixXd::Constant(M + 2, M + 2, 1);
  v.col(0) = Eigen::VectorXd::Zero(M + 2);
  v.col(M + 1) = Eigen::VectorXd::Zero(M + 2);
  v.row(0) = Eigen::VectorXd::Zero(M + 2);
  v.row(M + 1) = Eigen::VectorXd::Zero(M + 2);
  for (unsigned int i = 1; i < M + 1; ++i) {
    for (unsigned int j = 1; j < M + 1; ++j) {
      if (i <= midpoint && j >= midpoint) {
        v(i, j) = 0;
      }
    }
  }
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
  std::vector<unsigned> L_vec({6, 7, 8, 9,10,11});
  // ************************************************************
  // Code to be supplemented
  // ************************************************************
}
/* SAM_LISTING_END_9 */

}  // namespace MFMGLshape
