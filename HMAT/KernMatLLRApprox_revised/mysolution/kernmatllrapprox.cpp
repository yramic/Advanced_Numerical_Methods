/**
 * @file kernmatllrapprox.cpp
 * @brief NPDE homework KernMatLLRApprox code
 * @author R. Hiptmair
 * @date September 2023
 * @copyright Developed at SAM, ETH Zurich
 */

#include "kernmatllrapprox.h"

#include <Eigen/src/Core/Matrix.h>

#include <chrono>
#include <cstddef>
#include <limits>

namespace KernMatLLRApprox {

/* SAM_LISTING_BEGIN_1 */
bool validateLLR(unsigned int q, double tol, double eta) {
  const int npts = 64;  // Number of collocation points
  // **********************************************************************
  // TO BE SUPPLEMENTED
  // **********************************************************************
  return true;
}

/* SAM_LISTING_END_1 */

// 1D Logarithmic kernel function
/* SAM_LISTING_BEGIN_2 */
struct LogKernel {
  LogKernel() = default;
  double operator()(double x, double y) const {
    const double d = std::abs(x - y);
    return ((d > (std::numeric_limits<double>::epsilon() * std::abs(x + y)))
                ? (-std::log(d))
                : 0.0);
  }
};

void tabulateConvergenceLLR(std::vector<unsigned int> &&n_vec,
                            std::vector<unsigned int> &&q_vec, double eta) {
  LogKernel G;
  Eigen::MatrixXd relerr(n_vec.size(), q_vec.size());
  for (int n_idx = 0; n_idx < n_vec.size(); ++n_idx) {
    for (int q_idx = 0; q_idx < q_vec.size(); ++q_idx) {
      // Create equidistant points
      std::vector<HMAT::Point<1>> pts;
      for (int pt_idx = 0; pt_idx < n_vec[n_idx]; pt_idx++) {
        HMAT::Point<1> p;
        p.idx = pt_idx;
        p.x[0] = static_cast<double>(pt_idx) / (n_vec[n_idx] - 1);
        pts.push_back(p);
      }
      // **********************************************************************
      // TO BE SUPPLEMENTED
      // **********************************************************************
    }
  }

  std::cout << "Relative errors of local low-rank matrix approximation"
            << std::endl;
  std::cout << "n\\q ";
  for (int q_idx = 0; q_idx < q_vec.size(); ++q_idx) {
    std::cout << std::setw(15) << q_vec[q_idx];
  }
  std::cout << std::endl;
  for (int n_idx = 0; n_idx < n_vec.size(); ++n_idx) {
    std::cout << std::setw(15) << n_vec[n_idx];
    for (int q_idx = 0; q_idx < q_vec.size(); ++q_idx) {
      std::cout << std::setw(15) << relerr(n_idx, q_idx);
    }
    std::cout << std::endl;
  }
}
/* SAM_LISTING_END_2 */

/* SAM_LISTING_BEGIN_3 */
void runtimeMatVec(std::vector<unsigned int> &&n_vec, unsigned int n_runs,
                   unsigned int q, double eta) {
  std::cout << "runtimeMatVec(n_runs = " << n_runs << ", q = " << q
            << ", eta = " << eta << ")" << std::endl;
  LogKernel G;
  unsigned int n_cnt = n_vec.size();
  for (unsigned int n : n_vec) {
    // Create equidistant points
    std::vector<HMAT::Point<1>> pts;
    for (int pt_idx = 0; pt_idx < n; pt_idx++) {
      HMAT::Point<1> p;
      p.idx = pt_idx;
      p.x[0] = static_cast<double>(pt_idx) / (n - 1);
      pts.push_back(p);
    }
    // Allocate cluster tree objects (the same for both directions)
    auto T_row = std::make_shared<KernMatLLRApprox::LLRClusterTree>(q, pts);
    auto T_col = std::make_shared<KernMatLLRApprox::LLRClusterTree>(q, pts);
    // Build local low-rank compressed matrix
    KernMatLLRApprox::BiDirChebPartMat1D<LogKernel> Mt(T_row, T_col, G, q, eta);
// **********************************************************************
// TO BE SUPPLEMENTED
// **********************************************************************
  }
}
/* SAM_LISTING_END_3 */

}  // namespace KernMatLLRApprox
