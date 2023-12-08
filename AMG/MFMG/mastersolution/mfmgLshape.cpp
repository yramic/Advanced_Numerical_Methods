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
  double midpoint = (M + 1) / 2;
  for (unsigned int i = 1; i < M + 1; ++i) {
    for (unsigned int j = 1; j < M + 1; ++j) {
      //Evaluating the Matrix-vector product directly
      eta(i, j) = (4.0 + h*h * c_) * mu(i, j) +
                  (-mu(i - 1, j) - mu(i + 1, j) - mu(i, j - 1) - mu(i, j + 1));
      if (i <= midpoint && j >= midpoint) {
        eta(i, j) = 0;
      }
    }
  }
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

  // We need 4 (spatial) indices: two identify the hat function in space that is in the domain of $\cob{A}$,
  // two different ones identify the hat function in the range of $\cob{A}$.
  // (the dimension of $\cob{A}$ is $M*M\times M*M$, the following indices go from $0$ to $M-1$ respectively.)
  int ir = 0;
  int jr = 0;
  int ic = 0;
  int jc = 0;
  for (int k = 0; k < A.outerSize(); ++k) {
    for (Eigen::SparseMatrix<double>::InnerIterator it(A, k);
         it; ++it) {
      ir = it.row() % M;  // Row index of the unit square
      jr = it.row() / M;  // Column index of the unit square
      ic = it.col() % M;  // Row index of the unit square
      jc = it.col() / M;  // Column index of the unit square
      if ((ir <= midpoint - 1 && jr >= midpoint - 1) ||
          (ic <= midpoint - 1 && jc >= midpoint - 1)) {
        // We are in the upper right quadrant
        if (it.col() == it.row()) {
          // Diagonal entry
          it.valueRef() = 1.;
        }
        if (it.col() != it.row()) {
          // Entries away from the diagonal
          it.valueRef() = 0.;
        }
      }
    }
  }
  Eigen::SparseLU<Eigen::SparseMatrix<double>> sparse_solver;
  sparse_solver.analyzePattern(A);
  sparse_solver.factorize(A);
  assert(sparse_solver.info() == Eigen::Success);
  sol = sparse_solver.solve(phi_v);
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
  for (unsigned int i = 1; i < M + 1; ++i) {
    for (unsigned int j = 1; j < M + 1; ++j) {
      mu(i, j) = (phi(i, j) +
                  (mu(i - 1, j) + mu(i + 1, j) + mu(i, j - 1) + mu(i, j + 1))) /
                 (4.0 + h * h * c_);
      if (i <= midpoint && j >= midpoint) {
        mu(i, j) = 0;
      }
    }
  }
}
/* SAM_LISTING_END_3 */

/* SAM_LISTING_BEGIN_4 */
GridFunction DirichletBVPMultiGridSolver::residual(
    unsigned int level, const GridFunction &mu, const GridFunction &phi) const {
  const unsigned int M = std::pow(2, level) - 1;
  GridFunction res(M + 2, M + 2);
  res = (phi - applyGridOperator(level, mu));
  return res;
}

/* SAM_LISTING_END_4 */

/* SAM_LISTING_BEGIN_5 */
GridFunction DirichletBVPMultiGridSolver::prolongate(
    unsigned int level, const GridFunction &gamma) const {
  const unsigned int M = std::pow(2, level) - 1;
  GridFunction zeta((M + 2), (M + 2));
  zeta.setZero();
  for (unsigned int i = 1; i < M + 1; ++i) {
    //odd
    for (unsigned int j = 1; j < M + 1; ++j) {
      //odd j
      zeta(i, j) =
          0.25 *
          (gamma((i - 1) / 2, (j - 1) / 2) + gamma((i - 1) / 2, (j + 1) / 2) +
           gamma((i + 1) / 2, (j - 1) / 2) + gamma((i + 1) / 2, (j + 1) / 2));
      // even j
      ++j;
      zeta(i, j) =
          0.5 * (gamma((i - 1) / 2, j / 2) + gamma((i + 1) / 2, j / 2));
    }
    ++i;
    // even i
    for (unsigned int j = 1; j < M + 1; ++j) {
      //odd j
      zeta(i, j) =
          0.5 * (gamma(i / 2, (j - 1) / 2) + gamma(i / 2, (j + 1) / 2));
      // even j
      ++j;
      zeta(i, j) = gamma(i / 2, j / 2);
    }
  }
  return zeta;
}

/* SAM_LISTING_END_5 */

/* SAM_LISTING_BEGIN_6 */
GridFunction DirichletBVPMultiGridSolver::restrict(
    unsigned int level, const GridFunction &rho) const {
  const unsigned int M = std::pow(2, level - 1) - 1;
  GridFunction sigma((M + 2), (M + 2));
  sigma.setZero();
  for (unsigned int i = 1; i < M + 1; ++i) {
    for (unsigned int j = 1; j < M + 1; ++j) {
      sigma(i, j) = rho(2 * i, 2 * j);
      sigma(i, j) += 0.5 * (rho(2 * i, 2 * j + 1) + rho(2 * i, 2 * j - 1) +
                            rho(2 * i + 1, 2 * j) + rho(2 * i - 1, 2 * j));
      sigma(i, j) +=
          0.25 * (rho(2 * i - 1, 2 * j - 1) + rho(2 * i - 1, 2 * j + 1) +
                  rho(2 * i + 1, 2 * j - 1) + rho(2 * i + 1, 2 * j + 1));
    }
  }
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

  unsigned counter = 0;
  // applies pre-smoothers and restrictions until L0 is reached
  for (unsigned l = L_; l > L0; --l) {
    sweepGaussSeidel(l, mu_vec.back(), phi_vec.back());  // Pre smoothing
    residual_vec.push_back(
        residual(l, mu_vec.back(),
                 phi_vec.back()));  // Compute residual and adds to residual
    phi_vec.push_back(restrict(
        l, residual_vec.back()));  // Restrict residual and append 
    GridFunction zero = phi_vec.back();
    zero.setZero();
    mu_vec.push_back(zero);  // Adds the zero grid function
    ++counter;
  }
  mu_vec.back() = directSolve(L0, phi_vec.back());
  for (unsigned l = L0; l < L_; ++l) {
    mu_vec[counter - 1] += prolongate(l + 1, mu_vec[counter]);  //prolongate mu
    sweepGaussSeidel(l + 1, mu_vec[counter - 1],
                     phi_vec[counter - 1]);  //post-smoothening
    --counter;
  }
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
  // Power iteration
  int maxiter = 10;
  int iter = 0;
  do {
    lambda_old = lambda_new;
    v /= v.norm();
    v_old = v;
    Av = solver.applyGridOperator(L, v);
    v = solver.multigridIteration(0 * v, Av, 2);
    v = v_old - v;
    lambda_new = v.norm();
    iter++;

  } while (std::abs(lambda_new - lambda_old) > tol * lambda_new &&
           iter < maxiter);
  return lambda_new;
}
/* SAM_LISTING_END_8 */

/* SAM_LISTING_BEGIN_9 */
void tabulateMGConvergenceRate() {
  std::vector<double> c_vec({-40, -20, -10, -1, 0, 1, 10, 20, 40});
  std::vector<unsigned> L_vec({6, 7, 8, 9,10,11});
  std::cout << "c/L";
  for (auto L : L_vec) {
    std::cout << std::setw(15) << L;
  }
  std::cout << std::endl;
  for (auto c : c_vec) {
    std::cout << c;
    for (auto L : L_vec) {
      double lambda = estimateMGConvergenceRate(c, L, 0.01);
      std::cout << std::setw(15) << lambda;
    }
    std::cout << std::endl;
  }
}
/* SAM_LISTING_END_9 */

}  // namespace MFMGLshape
