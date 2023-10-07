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
#if SOLUTION
  // Initialize collocation points (the same in x/y-direction)
  std::vector<HMAT::Point<1>> pts;
  for (int n = 0; n < npts; n++) {
    HMAT::Point<1> p;
    p.idx = n;
    p.x[0] = static_cast<double>(n) / (npts - 1);
    pts.push_back(p);
  }
  // Allocate cluster tree object (the same for both directions)
  auto T_row =
      std::make_shared<KernMatLLRApprox::LLRClusterTree<HMAT::CtNode<1>>>(q);
  T_row->init(pts);
  auto T_col =
      std::make_shared<KernMatLLRApprox::LLRClusterTree<HMAT::CtNode<1>>>(q);
  T_col->init(pts);

  // Loop over polynomials up to degree q and check whether for a polynomial
  // kernel the approximation error will vanish.
  for (unsigned int k = 0; k < q; ++k) {
    for (unsigned int l = 0; l < q; ++l) {
      struct PolynomialKernel {
        explicit PolynomialKernel() : deg_x_(0), deg_y_(0) {}
        PolynomialKernel(unsigned int deg_x, unsigned int deg_y)
            : deg_x_(deg_x), deg_y_(deg_y) {}
        double operator()(double x, double y) const {
          return std::pow(x, deg_x_) * std::pow(y, deg_y_);
        }
        unsigned int deg_x_;
        unsigned int deg_y_;
      } G(k, l);
      // Initialize local low-rank compressed matrix data structure
      //
      KernMatLLRApprox::BiDirChebPartMat1D<PolynomialKernel> Mt(T_row, T_col, G,
                                                                q, eta);
      std::cout << "validateLLR: BlockPartition: " << Mt.ffb_cnt
                << " far field blocks, " << Mt.nfb_cnt << " near-field blocks"
                << std::endl;
      auto [error, norm] = approxErrorLLR<PolynomialKernel>(Mt);
      std::cout << "Kernel (x,y) -> x^" << k << " * y^" << l
                << ": rel error = " << error / norm << std::endl;
      if (error > tol * norm) {
        return false;
      }
    }
  }
  return true;
#else
  // **********************************************************************
  // TO BE SUPPLEMENTED
  // **********************************************************************
  return true;
#endif
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
#if SOLUTION
      // Allocate cluster tree objects (the same for both directions)
      auto T_row =
          std::make_shared<KernMatLLRApprox::LLRClusterTree<HMAT::CtNode<1>>>(
              q_vec[q_idx]);
      T_row->init(pts);
      auto T_col =
          std::make_shared<KernMatLLRApprox::LLRClusterTree<HMAT::CtNode<1>>>(
              q_vec[q_idx]);
      T_col->init(pts);
      // Build local low-rank compressed matrix
      KernMatLLRApprox::BiDirChebPartMat1D<LogKernel> Mt(T_row, T_col, G,
                                                         q_vec[q_idx], eta);
      std::cout << "tabulateConvergenceLLR: q=" << q_vec[q_idx]
                << ", n_pts = " << n_vec[n_idx]
                << " : BlockPartition: " << Mt.ffb_cnt << " far field blocks, "
                << Mt.nfb_cnt << " near-field blocks" << std::endl;
      auto [error, norm] = approxErrorLLR<LogKernel>(Mt);
      std::cout << "abs err = " << error << ", rel err = " << error / norm
                << std::endl;
      relerr(n_idx, q_idx) = error / norm;
#else
      // **********************************************************************
      // TO BE SUPPLEMENTED
      // **********************************************************************
#endif
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
    auto T_row =
        std::make_shared<KernMatLLRApprox::LLRClusterTree<HMAT::CtNode<1>>>(q);
    T_row->init(pts);
    auto T_col =
        std::make_shared<KernMatLLRApprox::LLRClusterTree<HMAT::CtNode<1>>>(q);
    T_col->init(pts);
    // Build local low-rank compressed matrix
    KernMatLLRApprox::BiDirChebPartMat1D<LogKernel> Mt(T_row, T_col, G, q, eta);
#if SOLUTION
    const size_t nrows = Mt.rows();
    const size_t ncols = Mt.cols();
    Eigen::VectorXd x = Eigen::VectorXd::Constant(ncols, 1.0);
    Eigen::VectorXd y(nrows);
    // Average runtimes over n\_run runs
    double ms_time = std::numeric_limits<double>::max();
    for (int r = 0; r < n_runs; ++r) {
      auto t1 = std::chrono::high_resolution_clock::now();
      y = mvLLRPartMat(Mt, x);
      auto t2 = std::chrono::high_resolution_clock::now();
      /* Getting number of milliseconds as a double. */
      std::chrono::duration<double, std::milli> ms_double = (t2 - t1);
      ms_time = std::min(ms_time, ms_double.count());
    }
    std::cout << "n = " << n << ": " << ms_time << " ms, #nfb = " << Mt.nfb_cnt
              << ", #ffb = " << Mt.ffb_cnt << std::endl;
#else
// **********************************************************************
// TO BE SUPPLEMENTED
// **********************************************************************
#endif
  }
}
/* SAM_LISTING_END_3 */

}  // namespace KernMatLLRApprox
