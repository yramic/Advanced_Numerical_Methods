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
void outIPNode(const InterpNode<1> &ipnode, bool printV, std::ostream &o) {
  o << "IPNode: ";
  for (const HMAT::Point<1> pt : ipnode.pts) {
    o << "(" << pt.idx << ", " << pt.x[0] << ") ";
  }
  if (printV) {
    o << ", V = " << std::endl << ipnode.V;
  }
  o << std::endl;
}

/* SAM_LISTING_BEGIN_1 */


// One way to compute the error is with a FUNCTOR another would be to use  
// a lambda function
class KernelStrct {
  private:
    int R, C;    
  public:
    KernelStrct() = default; // Default constructor
    KernelStrct(int R_, int C_) : R(R_), C(C_) {} // Constructor
    // Functor (Function Object)
    double operator()(double x, double y) const {
      return std::pow(x,R)*std::pow(y,C);
    }  
};

bool validateLLR(unsigned int q, double tol, double eta) {
  const int npts = 64;  // Number of collocation points
  // **********************************************************************
  // TO BE SUPPLEMENTED 2-4k:
  // The task is to check if the compression error for 64 test points is indeed
  // zero!

  // First we need to setup the collocation points (test points):
  std::vector<HMAT::Point<1>> pts; // Note: Points are in 1D again
  for (unsigned int i {0}; i < npts; ++i) {
    HMAT::Point<1> pt;
    pt.idx = i;
    // Note points need to be equidistant:
    // Since the result is a double an npts as well as i are integers, a static
    // cast is necessary!
    pt.x[0] = static_cast<double> (i) / (npts - 1);
    pts.push_back(pt);
  }

  // Next: Allocate cluster tree object (the same for both directions, due to
  // the fact that the collocation points are equidistant)
  auto T_row = std::make_shared<
      KernMatLLRApprox::LLRClusterTree<KernMatLLRApprox::InterpNode<1>>>(q);
  T_row->init(pts);

  auto T_col = std::make_shared<
      KernMatLLRApprox::LLRClusterTree<KernMatLLRApprox::InterpNode<1>>>(q);
  T_col->init(pts);

  // Loop over polynomials up to degree q and check whether for a polynomial
  // kernel the approximation error will vanish.

  for (unsigned int k {0}; k < q; k++) {
    for (unsigned int l {0}; l < q; l++) {
      KernelStrct G(k,l);
      // decltype():
      // Inspects the declared type of an entity or an expression. 
      KernMatLLRApprox::BiDirChebPartMat1D<decltype(G)> Mt(T_row, T_col, G, q, eta);
      auto errors = approxErrorLLR<decltype(G)>(Mt);
      // errors.first = error
      // error.second = norm
      // std::cout << "The error is: " << errors.first << std::endl;
      // std::cout << "The norm is: " << errors.second << std::endl;
      if (errors.first > tol * errors.second) return false;
    }
  }
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
      // TO BE SUPPLEMENTED PROBLEM 2-4L:
      // Relative Compression Error needs to be computed for an asymptotically
      // smooth logarithmic kernel function

      // The Problem here is similar to the task before:
      // Again: Allocate cluster tree object (the same for both directions)
      // Also here it is actually the same since the collocation points are equidistant:
      auto T_row = std::make_shared<KernMatLLRApprox::LLRClusterTree<KernMatLLRApprox::InterpNode<1>>>(q_vec[q_idx]);
      T_row->init(pts);

      auto T_col = std::make_shared<KernMatLLRApprox::LLRClusterTree<KernMatLLRApprox::InterpNode<1>>>(q_vec[q_idx]);
      T_col->init(pts);

      // Next, again: Build local the low-rank compressed matrix
      KernMatLLRApprox::BiDirChebPartMat1D<decltype(G)> Mt(T_row, T_col, G, q_vec[q_idx], eta);
      
      // Compute the error and norm again:
      auto errors = approxErrorLLR<decltype(G)>(Mt);
      // errors.first = error
      // errors.second = norm
      relerr(n_idx,q_idx) = errors.first/errors.second;

      // std::cout << errors.first << std::endl;

      // The error gets not printed, so far I couldn't find the error, when I run
      // the test and tried to run some more tests, it stated that Segmentation fault
      // So, setting up the cluster tree is not working properly
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
    auto T_row = std::make_shared<
        KernMatLLRApprox::LLRClusterTree<KernMatLLRApprox::InterpNode<1>>>(q);
    T_row->init(pts);
    auto T_col = std::make_shared<
        KernMatLLRApprox::LLRClusterTree<KernMatLLRApprox::InterpNode<1>>>(q);
    T_col->init(pts);
    // Build local low-rank compressed matrix
    KernMatLLRApprox::BiDirChebPartMat1D<LogKernel> Mt(T_row, T_col, G, q, eta);
// **********************************************************************
// TO BE SUPPLEMENTED PROBLEM 2-4M:
  const size_t nrows {Mt.rows()};
  const size_t ncols {Mt.cols()};

  Eigen::VectorXd x(ncols), y(nrows);
  x = Eigen::VectorXd::Ones(ncols);

  double duration {0};
  for (unsigned int i{0}; i < n_runs; ++i) {
    auto t1 = std::chrono::high_resolution_clock::now();
    y = mvLLRPartMat(Mt, x);
    auto t2 = std::chrono::high_resolution_clock::now();
    // Getting the number of miliseconds:
    std::chrono::duration<double,std::milli> ms_double = t2 - t1;
    duration = ( (duration == 0) ? ms_double.count() : std::min(duration, ms_double.count())); 
  }

  std::cout << "n: " << n << "; Runtime: " << duration << " ms" << std::endl;
// **********************************************************************
  }
}
/* SAM_LISTING_END_3 */

}  // namespace KernMatLLRApprox
