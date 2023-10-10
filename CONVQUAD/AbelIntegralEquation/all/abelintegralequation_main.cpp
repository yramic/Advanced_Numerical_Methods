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

#if SOLUTION
    cout << "\nSpectral Galerkin\n" << endl;
    double err_max, err_max_alt;
    for (int p = 2; p <= 10; ++p) {
      VectorXd u_app = AbelIntegralEquation::poly_spec_abel(y, p, tau);
      VectorXd diff = u_ex - u_app;
      err_max = diff.cwiseAbs().maxCoeff();
      double dp = p;
      if (p == 2) {
        cout << "p = " << p << setw(15) << "Max = " << scientific
             << setprecision(3) << err_max << endl;
      } else {
        cout << "p = " << p << setw(15) << "Max = " << scientific
             << setprecision(3) << err_max << setw(15) << " EOC = "
             << std::log2(err_max_alt / err_max) / std::log2((dp + 1) / dp)
             << endl;
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
      VectorXd grid = VectorXd::LinSpaced(N + 1, 0., 1.);
      VectorXd u_ex(N + 1);
      for (int i = 0; i < N + 1; ++i) {
        u_ex(i) = u(grid(i));
      }

      VectorXd u_app = AbelIntegralEquation::cq_ieul_abel(y, N);
      VectorXd diff = u_ex - u_app;
      err_max = diff.cwiseAbs().maxCoeff();
      if (N == 16) {
        cout << "N = " << N << setw(15) << "Max = " << scientific
             << setprecision(3) << err_max << endl;
      } else {
        cout << "N = " << N << setw(15) << "Max = " << scientific
             << setprecision(3) << err_max << setw(15)
             << " EOC = " << std::log2(err_max_alt / err_max) << endl;
      }

      err_max_alt = err_max;
    }

    cout << "\n\nConvolution Quadrature, BDF-2\n" << endl;
    for (int N = 16; N <= 2048; N *= 2) {
      VectorXd grid = VectorXd::LinSpaced(N + 1, 0., 1.);
      VectorXd u_ex(N + 1);
      for (int i = 0; i < N + 1; ++i) {
        u_ex(i) = 2. / M_PI * sqrt(grid(i));
      }

      VectorXd u_app = AbelIntegralEquation::cq_bdf2_abel(y, N);
      VectorXd diff = u_ex - u_app;
      err_max = diff.cwiseAbs().maxCoeff();
      if (N == 16) {
        cout << "N = " << N << setw(15) << "Max = " << scientific
             << setprecision(3) << err_max << endl;
      } else {
        cout << "N = " << N << setw(15) << "Max = " << scientific
             << setprecision(3) << err_max << setw(15)
             << " EOC = " << std::log2(err_max_alt / err_max) << endl;
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