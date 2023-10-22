/**
 * @file stationarylineariterations.cpp
 * @brief NPDE homework StationaryLinearIterations code
 * @author D. Casati
 * @date October 2018
 * @copyright Developed at SAM, ETH Zurich
 */

#include "stationarylineariterations.h"

namespace StationaryLinearIterations {
/* SAM_LISTING_BEGIN_1 */
using triplet = Eigen::Triplet<double>;
Eigen::SparseMatrix<double> poissonMatrix(unsigned int n) {
  int N = (n - 1) * (n - 1);  // Matrix size
  // define vector of triplets and reserve memory
  std::vector<triplet> entries;  // For temporary COO format
  // set the vector of triplets
  for (int block_id = 0; block_id < (n - 1); block_id++) {
    int start_id = block_id * (n - 1);
    int end_id = (block_id + 1) * (n - 1) - 1;
    // tri-diagonal matrix T
    for (int i = start_id; i <= end_id; i++) {
      entries.emplace_back(i, i, 4);  // Diagonal entry
      if (i > start_id) {             // T(i-1,i) = T(i,i-1)
        entries.emplace_back(i, i - 1, -1.0);
        entries.emplace_back(i - 1, i, -1.0);
      }
    }
    // identity blocks
    if (block_id > 0) {
      for (int i = start_id; i <= end_id; i++) {
        entries.emplace_back(i, i - (n - 1), -1.0);
        entries.emplace_back(i - (n - 1), i, -1.0);
      }
    }
  }
  // create the sparse matrix in CCS format
  Eigen::SparseMatrix<double> A(N, N);
  A.setFromTriplets(entries.begin(), entries.end());
  return A;
}
/* SAM_LISTING_END_1 */

/* SAM_LISTING_BEGIN_2 */
void gaussSeidel(const Eigen::SparseMatrix<double> &A,
                 const Eigen::VectorXd &phi, Eigen::VectorXd &mu, int maxItr,
                 double TOL) {
  Eigen::VectorXd delta(A.rows());
  do {
    // Rely on Eigen's ability to solve a sparse triangular system efficiently
    delta = A.triangularView<Eigen::Lower>().solve(
        phi - A * mu);  // Triangular solve \Label[line]{sligs:slv}
    mu += delta;
  } while (delta.norm() > TOL * mu.norm());
}
/* SAM_LISTING_END_2 */

/* SAM_LISTING_BEGIN_4 */
double comp_lmax_gaussSeidel(const Eigen::SparseMatrix<double> &X,
                             double TOL = 1.0E-03) {
  int N = X.rows();

  Eigen::VectorXd v = Eigen::VectorXd::Random(N);
  double lambda_new = 0;
  double lambda_old = 1;

  int itr = 0;
  // Power iteration 
  do {
    lambda_old = lambda_new;
    v /= v.norm();
    v -= X.triangularView<Eigen::Lower>().solve(X * v);
    lambda_new = v.norm();
  } while (fabs(lambda_new - lambda_old) > TOL * lambda_new);

  return lambda_new;
}
/* SAM_LISTING_END_4 */

/* SAM_LISTING_BEGIN_3 */
double gaussSeidelRate(unsigned int n, double c, double TOL) {
  int N = (n - 1) * (n - 1);

  Eigen::SparseMatrix<double> X = poissonMatrix(n);
  Eigen::SparseMatrix<double> I(N, N);
  I.setIdentity();
  // Is this efficient?
  X += c * I;
  double lambda_max = comp_lmax_gaussSeidel(X, TOL);
  return lambda_max;
}
/* SAM_LISTING_END_3 */

}  // namespace StationaryLinearIterations
