#include "abelintegralequation.h"

int main() {
  /* SAM_LISTING_BEGIN_1 */
  {
    auto u = [](double t) { return 2. / M_PI * sqrt(t); };
    auto y = [](double t) { return t; };

    double tau = 0.01;
    const std::size_t N = round(1. / tau);
    const Eigen::VectorXd grid = Eigen::VectorXd::LinSpaced(N + 1, 0., 1.);
    Eigen::VectorXd u_ex(N + 1);
    for (int i = 0; i < N + 1; ++i) {
      u_ex(i) = u(grid(i));
    }

#if SOLUTION
    std::cout << "\nSpectral Galerkin\n" << std::endl;
    double err_max, err_max_alt;
    for (int p = 2; p <= 10; ++p) {
      Eigen::VectorXd u_app = AbelIntegralEquation::poly_spec_abel(y, p, tau);
      Eigen::VectorXd diff = u_ex - u_app;
      err_max = diff.cwiseAbs().maxCoeff();
      double dp = p;
      if (p == 2) {
        std::cout << "p = " << p << std::setw(15) << "Max = " << std::scientific
             << std::setprecision(3) << err_max << std::endl;
      } else {
        std::cout << "p = " << p << std::setw(15) << "Max = " << std::scientific
             << std::setprecision(3) << err_max << std::setw(15) << " EOC = "
             << std::log2(err_max_alt / err_max) / std::log2((dp + 1) / dp)
             << std::endl;
      }
      err_max_alt = err_max;
    }
#else
// **********************************************************************
// Your Solution here
// **********************************************************************/
#endif
  }
  /* SAM_LISTING_END_1 */

  /* SAM_LISTING_BEGIN_4 */

#if SOLUTION
  {
    auto u = [](double t) { return 2. / M_PI * sqrt(t); };
    auto y = [](double t) { return t; };

    double err_max, err_max_alt;
    cout << "\n\nConvolution Quadrature, Implicit Euler\n" << endl;
    for (int N = 16; N <= 2048; N *= 2) {
      Eigen::VectorXd grid = Eigen::VectorXd::LinSpaced(N + 1, 0., 1.);
      Eigen::VectorXd u_ex(N + 1);
      for (int i = 0; i < N + 1; ++i) {
        u_ex(i) = u(grid(i));
      }

      Eigen::VectorXd u_app = AbelIntegralEquation::cq_ieul_abel(y, N);
      Eigen::VectorXd diff = u_ex - u_app;
      err_max = diff.cwiseAbs().maxCoeff();
      if (N == 16) {
        std::cout << "N = " << N << std::setw(15) << "Max = " << std::scientific
             << std::setprecision(3) << err_max << std::endl;
      } else {
        std::cout << "N = " << N << std::setw(15) << "Max = " << std::scientific
             << std::setprecision(3) << err_max << std::setw(15)
             << " EOC = " << std::log2(err_max_alt / err_max) << std::endl;
      }

      err_max_alt = err_max;
    }

    std::cout << "\n\nConvolution Quadrature, BDF-2\n" << std::endl;
    for (int N = 16; N <= 2048; N *= 2) {
      Eigen::VectorXd grid = Eigen::VectorXd::LinSpaced(N + 1, 0., 1.);
      Eigen::VectorXd u_ex(N + 1);
      for (int i = 0; i < N + 1; ++i) {
        u_ex(i) = 2. / M_PI * sqrt(grid(i));
      }

      Eigen::VectorXd u_app = AbelIntegralEquation::cq_bdf2_abel(y, N);
      Eigen::VectorXd diff = u_ex - u_app;
      err_max = diff.cwiseAbs().maxCoeff();
      if (N == 16) {
        std::cout << "N = " << N << std::setw(15) << "Max = " << std::scientific
             << std::setprecision(3) << err_max << std::endl;
      } else {
        cout << "N = " << N << std::setw(15) << "Max = " << std::scientific
             << std::setprecision(3) << err_max << std::setw(15)
             << " EOC = " << std::log2(err_max_alt / err_max) << std::endl;
      }

      err_max_alt = err_max;
    }
  }
#else
// **********************************************************************
// Your Solution here
// **********************************************************************/
#endif
  /* SAM_LISTING_END_4 */
  return 0;
}