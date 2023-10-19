#include "abelintegralequation.h"

int main() {
  /* SAM_LISTING_BEGIN_1 */
  {
    // Exact solution
    auto u = [](double t) { return 2. / M_PI * sqrt(t); };
    auto y = [](double t) { return t; };

    const double tau = 0.01;
    const std::size_t N = std::round(1. / tau);
    const Eigen::VectorXd grid = Eigen::VectorXd::LinSpaced(N + 1, 0., 1.);
    // Exact solution at grid points
    const Eigen::VectorXd u_ex = Eigen::VectorXd::NullaryExpr(
        N + 1, [&](Eigen::Index i) { return u(grid(i)); });

    // **********************************************************************
    // Your Solution here
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