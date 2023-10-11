#include <Eigen/Dense>
#include <array>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <istream>

#include "gauleg.hpp"

//-----------------------------------------------------------------------------
struct TriaPanel {
  // Array of 3d-vectors containing the triangle's vertices
  std::array<Eigen::Vector3d, 3> v;
  TriaPanel(const Eigen::Vector3d& a, const Eigen::Vector3d& b,
            const Eigen::Vector3d& c) {
    // Initialize array of vertices with the corresponding points
    v[0] = a;
    v[1] = b;
    v[2] = c;
  }

  Eigen::Vector3d getVertex(int i) const {
    assert(i >= 0 && i <= 2);
    return v[i];
  }
};

//-----------------------------------------------------------------------------
/* SAM_LISTING_BEGIN_1 */
bool ProjectOnTria(const TriaPanel& T, const Eigen::Vector3d& x,
                   Eigen::Vector3d& xp) {
  // TODO: Implement your code
  return true;
}
/* SAM_LISTING_END_1 */

//-----------------------------------------------------------------------------
/* SAM_LISTING_BEGIN_3 */
double integrateTiSing(const Eigen::Vector2d& b, const Eigen::Vector2d& c,
                       double zeta, int n = 6) {
  double res = 0;
  // TODO: Implement your code
  return res;
}
/* SAM_LISTING_END_3 */

//-----------------------------------------------------------------------------
/* SAM_LISTING_BEGIN_5 */
double integrateTSing(const TriaPanel& T, const Eigen::Vector3d& x, int n = 6) {
  double res = 0.;
  // TODO: Implement your code
  return res;
}
/* SAM_LISTING_END_5 */

//------------------------------------------------------------------------------
/* SAM_LISTING_BEGIN_6 */
void testIntegrateTSing(const int& n) {
  // Test on reference triangle for x on vertices and edges (for which we
  // computed the exact solution using Maple).
  Eigen::Vector3d ah, bh, ch;
  ah << 0., 0., 0.;
  bh << 1., 0., 0.;
  ch << 0., 1., 0.;
  TriaPanel T0(ah, bh, ch);

  std::cout << " Testing for n = " << n << std::endl;

  std::cout << " Error on a :"
            << fabs(integrateTSing(T0, ah, n) + sqrt(2) * log(sqrt(2) - 1))
            << std::endl;

  std::cout << " Error on b :"
            << fabs(integrateTSing(T0, bh, n) + log(sqrt(2) - 1)) << std::endl;

  std::cout << " Error on c :"
            << fabs(integrateTSing(T0, ch, n) + log(sqrt(2) - 1)) << std::endl;

  std::cout << " Error on 0.5*(a+b) :";
  double exres = 1.676348272;
  std::cout << fabs(integrateTSing(T0, 0.5 * (ah + bh), n) - exres)
            << std::endl;

  std::cout << " Error on 0.5*(b+c) :";
  double exres2 = 0.5 * (log(sqrt(2) + 1) - 3 * log(sqrt(2) - 1));
  std::cout << fabs(integrateTSing(T0, 0.5 * (bh + ch), n) - exres2)
            << std::endl;

  std::cout << " Error on 0.5*(a+c) :";
  std::cout << fabs(integrateTSing(T0, 0.5 * (ah + ch), n) - exres) << std::endl
            << std::endl;
}
/* SAM_LISTING_END_6 */

//-----------------------------------------------------------------------------
int main() {
  // 1. Test convergence for integration over Ti
  Eigen::Vector2d b0, c0;
  b0 << 1., 0.;
  c0 << 0.5, 1.;
  double Iref = integrateTiSing(b0, c0, 0., 1000);
  Eigen::VectorXd error(19);
  Eigen::VectorXi N = Eigen::VectorXi::LinSpaced(19, 1, 19);
  // TODO: Implement your code

  // Output for plot
  std::ofstream out_errorQ("integrateTiSing_errors.txt");
  out_errorQ << std::setprecision(18) << error;
  out_errorQ.close();
  std::ofstream out_N("integrateTiSing_N.txt");
  out_N << N;
  out_N.close();

  // You may use testIntegrateTSing(n); to test integrateTSing for a toy case

  // 3. Test for given triangle moving the point x
  Eigen::Vector3d a, b, c;
  a << 1., 1., 1.;
  b << 2., 1., 0.;
  c << 0., 1., 0.;
  TriaPanel T(a, b, c);
  // TODO: Implement your code

  return 0;
}
