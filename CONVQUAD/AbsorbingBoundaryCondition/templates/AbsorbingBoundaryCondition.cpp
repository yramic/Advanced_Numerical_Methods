#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <Eigen/SparseCholesky>
#include <cmath>
#include <iomanip>
#include <iostream>

using namespace Eigen;
using namespace std;


VectorXd pconv(const VectorXd& u, const VectorXd& x) {
  using idx_t = VectorXd::Index; // may be unsigned !
  const idx_t n = x.size();
  VectorXd z = VectorXd::Zero(n);
  // Need signed indices when differences are formed
  for (long k = 0; k < n; ++k) {
      for (long j = 0; j < n; ++j) {
          long ind = (k - j < 0 ? n + k - j : k - j);
          z(k) += u(ind)*x(j);
      }
  }
  return z;
}


/* @brief Find the unknown function u at final time t = 1 in the evolution problem
 * using Galerkin discretization and convolution quadrature (BDF-2)
 * \param phi Template function for the right-hand side
 * \param M Number of discretization intervals in time
 * \param N Number of discretization intervals in space
 * \param p Order of quadrature rule
 * \\return Values of u at final time t = 1
 */
/* SAM_LISTING_BEGIN_0 */
template<typename FUNC>
VectorXd solveABC(const FUNC& phi, size_t M, size_t N, int p)
{
    // TODO: Find the unknown function u at final time t = 1 in the evolution problem
}
/* SAM_LISTING_END_0 */


int main() {
    /* SAM_LISTING_BEGIN_1 */
    // TODO: Tabulate the H1-error of the Galerkin discretization + convolution quadrature
    /* SAM_LISTING_END_1 */
}
