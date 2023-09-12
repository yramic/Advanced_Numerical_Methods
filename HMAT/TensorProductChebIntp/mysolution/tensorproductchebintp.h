/**
 * @file tensorproductchebintp.h
 * @brief NPDE homework TensorProductChebIntp code
 * @author R. Hiptmair
 * @date August 2023
 * @copyright Developed at SAM, ETH Zurich
 */

#include <Eigen/Dense>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <vector>

namespace TensorProductChebIntp {

/** @brief evaluation of Chybchev interpolant in many points
 * @tparam FUNCTOR a type compatible with std::function<double(double)>
 * @param q number of interpolation nodes = degree + 1
 * @param y sequence of data values for interpolation conditions at Chebychev
 * nodes
 * @param x vector of evaluation nodes
 *
 * Note that the number of y value agrees with $\cob{q}$, the polynomial degree
 * +1.
 */
/* SAM_LISTING_BEGIN_1 */
template <typename FUNCTOR>
std::vector<double> chebInterpEval1D(unsigned int q, FUNCTOR f,
                                     std::vector<double> &x) {
  // Initialize Chebychev nodes, compute barycentric weights
  // $\cob{\lambda_i}$ (up to q-dependent scaling) and sample the function
  std::vector<double> lambda(q);  // barycentric weights
  std::vector<double> t(q);       // Chebychev nodes
  std::vector<double> y(q);       // Sampled function values
  int sgn = 1;
  for (int k = 0; k < q; ++k, sgn *= -1) {
    t[k] = std::cos((2.0 * k - 1.0) / (2 * q) * M_PI);
    lambda[k] = sgn * std::sin((2.0 * k - 1.0) / (2 * q) * M_PI);
    y[k] = f(t[k]);
  }
  // Loop over all evaluation points
  std::vector<double>::size_type N = x.size();
  std::vector<double> res(N, 0.0);  // Result vector
  for (int k = 0; k < N; ++k) {
    // Denominator in the barycentric interpolation formula
    double den = 0.0;
    bool nonode = true;
    for (int j = 0; j < q; ++j) {
      if (x[k] == t[j]) {
        // Avoid division by zero
        res[k] = y[j];
        nonode = false;
        break;
      }
      const double tmp = lambda[j] / (x[k] - t[j]);
      res[k] += y[j] * tmp;
      den += tmp;
    }
    if (nonode) {
      res[k] /= den;
    }
  }
  return res;
}
/* SAM_LISTING_END_1 */

/* SAM_LISTING_BEGIN_2 */
template <typename FUNCTOR>
std::vector<double> chebInterpEval2D(unsigned int q, FUNCTOR f,
                                     std::vector<Eigen::Vector2d> &x) {
  // Initialize Chebychev nodes, compute barycentric weights
  // $\cob{\lambda_i}$ (up to q-dependent scaling) and sample the function
  std::vector<double> lambda(q);  // barycentric weights
  std::vector<double> t(q);       // Chebychev nodes
  int sgn = 1;
  for (int k = 0; k < q; ++k, sgn *= -1) {
    t[k] = std::cos((2.0 * k - 1.0) / (2 * q) * M_PI);
    lambda[k] = sgn * std::sin((2.0 * k - 1.0) / (2 * q) * M_PI);
  }
  Eigen::Matrix2d y(q, q);  // Sampled function values
  for (int k = 0; k < q; ++k) {
    for (int l = 0; l <= k; ++l) {
      y(k, l) = y(k, l) = f(t[k], t[l]);
    }
  }
  // Loop over all evaluation points
  std::vector<double>::size_type N = x.size();
  std::vector<double> res(N, 0.0);  // Result vector
  std::vector<double> wx(q);        // Weights $\cob{w_x^i}$
  std::vector<double> wy(q);        // Weights $\cob{w_y^i}$
  for (int k = 0; k < N; ++k) {
    // To be supplemented
  }
  return res;
}
/* SAM_LISTING_END_2 */

template <typename FUNCTOR>
std::vector<double> genChebInterpEval2D(unsigned int q, FUNCTOR f,
                                        Eigen::Vector2d a, Eigen::Vector2d b,
                                        std::vector<Eigen::Vector2d> &x) {}

}  // namespace TensorProductChebIntp
