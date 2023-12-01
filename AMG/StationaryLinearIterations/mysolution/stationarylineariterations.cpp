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
  // **********************************************************************
  // Problem 4-1:
  int idx_start, idx_end; // Variables to keep track of the start and end
  for (int idx_block{0}; idx_block < (n - 1); ++idx_block) {
    // This iterator goes through all "boxes" of the matrix A
    // Now we need to define the start and End for every Block:
    idx_start = idx_block * (n - 1);
    idx_end = (idx_block + 1) * (n - 1) - 1;
    // Tri-diagonal Matrix T:
    for (int i{idx_start}; i <= idx_end; ++i) {
      // First the Diagonal Entries:
      entries.emplace_back(i, i, 4.);
      if (i > idx_start) {
        // If this is the case then set the value above and on the left to -1!
        // In other words:
        // T(i-1, i) = T(i, i-1) = -1
        entries.emplace_back(i-1, i, -1.);
        entries.emplace_back(i, i-1, -1.);
      }
    }
    // For the Identity Block we actually use the same approach as before!
    // if idx_block > 0:
    // A(idx_block - 1, idx_block) = A(idx_block, idx_block -1) = -1
    if (idx_block > 0) {
      for (int i{idx_start}; i <= idx_end; ++i) {
        entries.emplace_back(i, i - (n-1), -1.);
        entries.emplace_back(i - (n-1), i, -1.);
      }
    } 
  }
  // **********************************************************************
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
  // **********************************************************************
  // Problem 4-1d:
  // Initialize the delta norm:
  Eigen::VectorXd delta(A.rows());
  // Counter for the number of iterations:
  unsigned int iteration{0};
  do {
    // Eigen can solve triangular systems efficiently!
    delta = A.triangularView<Eigen::Lower>().solve(phi - A * mu);
    mu += delta;
    ++iteration;
  } while (delta.norm() > TOL * mu.norm());
  
  // Print if convergence wasn't reached before iter == maxItr:
  if (iteration == maxItr) {
    std::cerr << "Did not converge to TOL: " << TOL
              << " in maximal number of iterations: " << iteration << std::endl;
  }
  // **********************************************************************
}
/* SAM_LISTING_END_2 */

/* SAM_LISTING_BEGIN_4 */
double comp_lmax_gaussSeidel(const Eigen::SparseMatrix<double> &X,
                             double TOL = 1.0E-03) {
  const int N = X.rows();

  Eigen::VectorXd v = Eigen::VectorXd::Random(N);
  double lambda_new = 0;
  double lambda_old = 1;

    // **********************************************************************
    // Problem 4-1e:
    // Here we need to find the max Eigenvalue of the Matrix X by using the 
    // Power Iteration Method from NumCSE 9.3.1!
    Eigen::VectorXd w(v.size());
    do {
      lambda_old = lambda_new;
      // It is worth noting that in the case of k = 0 the initial guess v is
      // just randomly assigned which can be seen above! Note that v is in my
      // notes the same as z!
      // Also it is based on the PSEUDOCODE 4.1.3.21
      v /= v.norm();
      v = X * v;
      lambda_new = v.norm();
    } while(std::abs(lambda_new - lambda_old) > TOL * lambda_new);
    // **********************************************************************
  return lambda_new;
}
/* SAM_LISTING_END_4 */

/* SAM_LISTING_BEGIN_3 */
double gaussSeidelRate(unsigned int n, double c, double TOL) {
  const unsigned int N = (n - 1) * (n - 1);

    // **********************************************************************
    // Problem 4-1e:
    Eigen::SparseMatrix<double> X = poissonMatrix(n);
    Eigen::SparseMatrix<double> I(N,N);
    I.setIdentity();
    X += c * I;
    const double lambda_max = comp_lmax_gaussSeidel(X, TOL);
    // **********************************************************************
    return lambda_max;
}
/* SAM_LISTING_END_3 */

}  // namespace StationaryLinearIterations
