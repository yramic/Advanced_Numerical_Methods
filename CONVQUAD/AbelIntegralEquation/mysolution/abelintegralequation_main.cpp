#include "abelintegralequation.h"

int main() {
  /* SAM_LISTING_BEGIN_1 */
  {
    // Exact solution
    auto u = [](double t) { return 2. / M_PI * sqrt(t); };
    auto y = [](double t) { return t; };

    // Generate points on the grid
    double tau = 0.01;
    size_t N = round(1. / tau);
    VectorXd grid = VectorXd::LinSpaced(N + 1, 0., 1.);
    // Exact solution at grid points
    VectorXd u_ex(N + 1);
    for (int i = 0; i < N + 1; ++i) {
      u_ex(i) = u(grid(i));
    }

// **********************************************************************
// PROBLEM 3-2g:
// Calculate the error between the approximated and the exact solution!
  std::cout << "\nSpectral Galerkin\n\n";
  double err_max, err_max_prior;
  for (int p{2}; p <= 10; ++p) {
    // Solution with Galerkin discretization with a polynomial basis
    const Eigen::VectorXd u_approx = AbelIntegralEquation::poly_spec_abel(y, p, tau);
    // Maximum Norm of discretization error:
    const Eigen::VectorXd diff = (u_ex - u_approx).cwiseAbs();
    // Compute error:
    // Note: for Eigen .begin() and .end() can't be used, instead .data() has to be used!
    err_max = *std::max_element(diff.data(), diff.data() + diff.size());

    std::cout << "p = " << p << std::setw(15) << "Max = " << std::scientific
              << std::setprecision(3) << err_max << std::setw(15);

    if (p > 2) {
      std::cout << " EOC = "
                << std::log2(err_max_prior / err_max) /
                    std::log2((static_cast<double>(p) + 1) / static_cast<double>(p));
    }
    std::cout << std::endl;
    err_max_prior = err_max;
  }
// **********************************************************************/
  }
  /* SAM_LISTING_END_1 */

  /* SAM_LISTING_BEGIN_4 */

// **********************************************************************
// Your Solution here
// **********************************************************************/
  /* SAM_LISTING_END_4 */
  return 0;
}