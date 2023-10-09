#include "abelintegralequation.h"

int main() {
  /* SAM_LISTING_BEGIN_1 */
  {
    auto u = [](double t) { return 2. / M_PI * sqrt(t); };
    auto y = [](double t) { return t; };

    double tau = 0.01;
    size_t N = round(1. / tau);
    VectorXd grid = VectorXd::LinSpaced(N + 1, 0., 1.);
    VectorXd u_ex(N + 1);
    for (int i = 0; i < N + 1; ++i) {
      u_ex(i) = u(grid(i));
    }

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