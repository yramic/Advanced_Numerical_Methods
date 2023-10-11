#include "abelintegralequation.h"

int main() {
  /* SAM_LISTING_BEGIN_1 */
  {
    auto u = [](double t) { return 2. / M_PI * sqrt(t); };
    auto y = [](double t) { return t; };

    const double tau = 0.01;
    const std::size_t N = std::round(1. / tau);
    const Eigen::VectorXd grid = Eigen::VectorXd::LinSpaced(N + 1, 0., 1.);
    const Eigen::VectorXd u_ex = Eigen::VectorXd::NullaryExpr(
        N + 1, [&](Eigen::Index i) { return u(grid(i)); });

    std::cout << "\nSpectral Galerkin\n\n";
    double err_max, err_max_alt;
    for (int p = 2; p <= 10; ++p) {
      const Eigen::VectorXd u_app =
          AbelIntegralEquation::poly_spec_abel(y, p, tau);
      const Eigen::VectorXd diff = u_ex - u_app;
      err_max = diff.cwiseAbs().maxCoeff();
      const double dp = p;

      std::cout << "p = " << p << std::setw(15) << "Max = " << std::scientific
                << std::setprecision(3) << err_max << std::setw(15);
      if (p > 2) {
        std::cout << " EOC = "
                  << std::log2(err_max_alt / err_max) /
                         std::log2((dp + 1) / dp);
      }
      std::cout << '\n';

      err_max_alt = err_max;
    }
  }
  /* SAM_LISTING_END_1 */

  /* SAM_LISTING_BEGIN_4 */

  {
    auto u = [](double t) { return 2. / M_PI * sqrt(t); };
    auto y = [](double t) { return t; };

    double err_max, err_max_alt;
    cout << "\n\nConvolution Quadrature, Implicit Euler\n\n";
    for (int N = 16; N <= 2048; N <<= 1) {
      const Eigen::VectorXd grid = Eigen::VectorXd::LinSpaced(N + 1, 0., 1.);
      const Eigen::VectorXd u_ex = Eigen::VectorXd::NullaryExpr(
          N + 1, [&](Eigen::Index i) { return u(grid(i)); });

      const Eigen::VectorXd u_app = AbelIntegralEquation::cq_ieul_abel(y, N);
      const Eigen::VectorXd diff = u_ex - u_app;
      err_max = diff.cwiseAbs().maxCoeff();

      std::cout << "N = " << N << std::setw(15) << "Max = " << std::scientific
                << std::setprecision(3) << err_max;
      if (N > 16) {
        std::cout << " EOC = " << std::log2(err_max_alt / err_max);
      }
      std::cout << '\n';

      err_max_alt = err_max;
    }

    std::cout << "\n\nConvolution Quadrature, BDF-2\n" << '\n';
    for (int N = 16; N <= 2048; N <<= 1) {
      Eigen::VectorXd grid = Eigen::VectorXd::LinSpaced(N + 1, 0., 1.);
      Eigen::VectorXd u_ex(N + 1);
      for (int i = 0; i < N + 1; ++i) {
        u_ex(i) = 2. / M_PI * sqrt(grid(i));
      }

      const Eigen::VectorXd u_app = AbelIntegralEquation::cq_bdf2_abel(y, N);
      const Eigen::VectorXd diff = u_ex - u_app;
      err_max = diff.cwiseAbs().maxCoeff();

      std::cout << "N = " << N << std::setw(15) << "Max = " << std::scientific
                << std::setprecision(3) << err_max;
      if (N > 16) {
        std::cout << " EOC = " << std::log2(err_max_alt / err_max);
      }
      std::cout << '\n';

      err_max_alt = err_max;
    }
  }
  /* SAM_LISTING_END_4 */
  return 0;
}