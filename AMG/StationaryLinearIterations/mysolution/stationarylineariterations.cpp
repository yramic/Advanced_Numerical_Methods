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
  // Code to be supplemented
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
  // Code to be supplemented
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
  // Code to be supplemented
  // **********************************************************************
  return lambda_new;
}
/* SAM_LISTING_END_4 */

/* SAM_LISTING_BEGIN_3 */
double gaussSeidelRate(unsigned int n, double c, double TOL) {
  const unsigned int N = (n - 1) * (n - 1);

  // **********************************************************************
  // Code to be supplemented
  // **********************************************************************
  return 0;
}
/* SAM_LISTING_END_3 */

}  // namespace StationaryLinearIterations
