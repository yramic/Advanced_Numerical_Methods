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
    Eigen::VectorXd grid = Eigen::VectorXd::LinSpaced(N + 1, 0., 1.);
    // Exact solution at grid points
    Eigen::VectorXd u_ex(N + 1);
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
// PROBLEM 3-2l:
{
    auto u = [](double t) { return 2. / M_PI * sqrt(t); };
    auto y = [](double t) { return t; };

    double err_max, err_max_alt;
    cout << "\n\nConvolution Quadrature, Implicit Euler\n" << endl;
  // The rhs function is defined by y(t) = t, for which the exact solution
  // can be computed as follows: u(t) = 2/pi * sqrt(t). These functions can
  // be expressed by lambda functions (Line 7 & 8)
  double err_cq_max, err_cq_max_prior;
  std::cout << "\n\nConvolution Quadrature, Implicit Euler\n\n";
  for (int N{16}; N <= 2048; N <<= 1) {
    // Note: This Loop doubles N each time!
    // Generate the Grid Points:
    const Eigen::VectorXd grid = Eigen::VectorXd::LinSpaced(N + 1, 0., 1.);
    // Exact solution at grid points
    Eigen::VectorXd u_ex(N + 1);
    for (int i = 0; i < N + 1; ++i) {
      u_ex(i) = u(grid(i));
    }
    // Solution for Convolution Quadrature based on implicit Euler method
    const Eigen::VectorXd u_approx = AbelIntegralEquation::cq_ieul_abel(y, N);
    // Get the Maximum Norm:
    const Eigen::VectorXd diff = (u_ex - u_approx).cwiseAbs(); // Absolute values
    err_cq_max = *std::max_element(diff.data(), diff.data() + diff.size()); // Find the max value
    std::cout << "N = " << N << std::setw(15) << "Max = " << std::scientific
              << std::setprecision(3) << err_cq_max << std::setw(15);
    
    if (N > 16) {
      std::cout << "EOC = " << std::log2(err_cq_max_prior / err_cq_max);
    }
    std::cout << std::endl;

    std::cout << "################## TEST ######################" << std::endl;
    std::cout << u_approx[0] << std::endl;
    std::cout << u_ex[0] << std::endl;
    std::cout << "################## TEST ######################" << std::endl;

    err_cq_max_prior = err_cq_max;
  }
  // Now same thing for the BDF-2 Method!
  std::cout << "\n\nConvolution Quadratre, BDF-2 Method\n\n";
  // Same Thing again:
  for (int N{16}; N <= 2048; N <<= 1) {
    // Generate Grid:
    const Eigen::VectorXd grid = Eigen::VectorXd::LinSpaced(N + 1, 0., 1.);
    // Exact solution at grid points
    Eigen::VectorXd u_ex(N + 1);
    for (int i = 0; i < N + 1; ++i) {
      u_ex(i) = u(grid(i));
    }
    // Solution using CQ based on BDF-2 Method
    const Eigen::VectorXd u_approx = AbelIntegralEquation::cq_bdf2_abel(y, N);
    // std::cout << u_approx[0] << std::endl;
    // std::cout << u_ex[0] << std::endl;
  //   // Maximum Norm of discretization error:
  //   const Eigen::VectorXd diff = (u_ex - u_approx).cwiseAbs();
  //   err_cq_max = *std::max(diff.data(), diff.data() + diff.size());
  //   std::cout << "N = " << N << std::setw(15) << "Max = " << std::scientific
  //             << std::setprecision(3) << err_cq_max << std::setw(15);
  //   if (N > 16) {
  //     std::cout << "EOC = " << std::log2(err_cq_max_prior / err_cq_max);
  //   }
  //   std::cout << std::endl;
  //   err_cq_max_prior = err_cq_max; // Update!
  }
}
  
// **********************************************************************/
  /* SAM_LISTING_END_4 */
  return 0;
}