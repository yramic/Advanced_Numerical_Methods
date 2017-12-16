#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <Eigen/SparseLU>
#include <cmath>
#include <iomanip>
#include <iostream>

using namespace Eigen;
using namespace std;


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
