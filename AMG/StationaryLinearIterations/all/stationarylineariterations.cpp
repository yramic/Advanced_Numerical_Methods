/**
 * @file stationarylineariterations.cpp
 * @brief ADVNCSE homework StationaryLinearIterations code
 * @author D. Casati, Bob Schreiner
 * @date October 2023
 * @copyright Developed at SAM, ETH Zurich
 */

#include "stationarylineariterations.h"

namespace StationaryLinearIterations {
/* SAM_LISTING_BEGIN_1 */
using triplet = Eigen::Triplet<double>;
Eigen::SparseMatrix<double> poissonMatrix(unsigned int n) {
  const int N = (n - 1) * (n - 1);  // Matrix size
  // define vector of triplets and reserve memory
  std::vector<triplet> entries;  // For temporary COO format
#if SOLUTION
  int start_id;
  int end_id;
  // set the vector of triplets

  for (int block_id = 0; block_id < (n - 1); block_id++) {
    start_id = block_id * (n - 1);
    end_id = (block_id + 1) * (n - 1) - 1;
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
#else
  // **********************************************************************
  // Code to be supplemented
  // **********************************************************************
#endif
  // create the sparse matrix in CRS format
  Eigen::SparseMatrix<double> A(N, N);
  A.setFromTriplets(entries.begin(), entries.end());
  A.makeCompressed();
  return A;
}
/* SAM_LISTING_END_1 */

/* SAM_LISTING_BEGIN_2 */
void gaussSeidel(const Eigen::SparseMatrix<double> &A,
                 const Eigen::VectorXd &phi, Eigen::VectorXd &mu, int maxItr,
                 double TOL) {
#if SOLUTION
  Eigen::VectorXd delta(A.rows());
  int iter = 0;
  do {
    // Rely on Eigen's ability to solve a sparse triangular system efficiently
    delta = A.triangularView<Eigen::Lower>().solve(
        phi - A * mu);  // Triangular solve \Label[line]{sligs:slv}
    mu += delta;
    ++iter;
  } while (delta.norm() > TOL * mu.norm() && iter < maxItr);
  if (maxItr == iter) {
    std::cerr << "Did not converge to tol: " << TOL
              << "in maximal number of iterations: " << iter << std::endl;
  }
#else
  // **********************************************************************
  // Code to be supplemented
  // **********************************************************************
#endif
}
/* SAM_LISTING_END_2 */

/* SAM_LISTING_BEGIN_4 */
double comp_lmax_gaussSeidel(const Eigen::SparseMatrix<double> &X,
                             double TOL = 1.0E-03) {
  const int N = X.rows();

  Eigen::VectorXd v = Eigen::VectorXd::Random(N);
  double lambda_new = 0;
  double lambda_old = 1;

#if SOLUTION
  // Power iteration
  do {
    lambda_old = lambda_new;
    v /= v.norm();
    v -= X.triangularView<Eigen::Lower>().solve(X * v);
    lambda_new = v.norm();
  } while (std::abs(lambda_new - lambda_old) > TOL * lambda_new);
#else
    // **********************************************************************
    // Code to be supplemented
    // **********************************************************************
#endif
  return lambda_new;
}
/* SAM_LISTING_END_4 */

/* SAM_LISTING_BEGIN_3 */
double gaussSeidelRate(unsigned int n, double c, double TOL) {
  const unsigned int N = (n - 1) * (n - 1);

#if SOLUTION
  Eigen::SparseMatrix<double> X = poissonMatrix(n);
  Eigen::SparseMatrix<double> I(N, N);
  I.setIdentity();
  X += c * I;
  const double lambda_max = comp_lmax_gaussSeidel(X, TOL);
  return lambda_max;
#else
    // **********************************************************************
    // Code to be supplemented
    // **********************************************************************
    return 0;
#endif
}
/* SAM_LISTING_END_3 */

}  // namespace StationaryLinearIterations
