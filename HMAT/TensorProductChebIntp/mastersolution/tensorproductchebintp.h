/**
 * @file tensorproductchebintp.h
 * @brief ADVNCSE homework TensorProductChebIntp code
 * @author R. Hiptmair , Bob Schreiner
 * @date August 2023
 * @copyright Developed at SAM, ETH Zurich
 */
#ifndef TENSORPRODCHEB_H_
#define TENSORPRODCHEB_H_

#include <Eigen/Dense>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <random>
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
  std::vector<double> lambda(q);  // barycentric weightsw
  std::vector<double> t(q);       // Chebychev nodes
  std::vector<double> y(q);       // Sampled function values
  int sgn = 1;
  for (int k = 0; k < q; ++k, sgn *= -1) {
    t[k] = std::cos((2.0 * k - 1.0) / (2 * q) * M_PI);  // \prbeqref{eq:chn}
    lambda[k] = sgn * std::sin((2.0 * k - 1.0) / (2 * q) * M_PI);  // \prbeqref{eq:bwf}
    y[k] = f(t[k]);                                        // $\cob{y_k}$
  }
  // Loop over all evaluation points $\cob{x_i}$
  const std::vector<double>::size_type N = x.size();
  std::vector<double> res(N, 0.0);  // Result vector
  for (int k = 0; k < N; ++k) {
    // Denominator in the barycentric interpolation formula
    double den = 0.0;
    bool nonode = true;
    for (int j = 0; j < q; ++j) {
      if (x[k] == t[j]) {
        // Avoid division by zero. In this case testing exact equality of
        // floating point numbers is safe, because the point of this test is not
        // to avoid amplification of roundoff errors, but really preventing a
        // division by zero exception.
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
                                     const std::vector<Eigen::Vector2d> &x) {
  // Initialize Chebychev nodes, compute barycentric weights
  // $\cob{\lambda_i}$ (up to q-dependent scaling) and sample the function
  std::vector<double> lambda(q);  // barycentric weights
  std::vector<double> t(q);       // Chebychev nodes
  int sgn = 1;
  for (int k = 0; k < q; ++k, sgn *= -1) {
    t[k] = std::cos((2.0 * k - 1.0) / (2 * q) * M_PI);
    lambda[k] = sgn * std::sin((2.0 * k - 1.0) / (2 * q) * M_PI);
  }
  Eigen::MatrixXd y(q, q);  // Sampled function values
  for (int k = 0; k < q; ++k) {
    for (int l = 0; l <= k; ++l) {
      y(k, l) = f(t[k], t[l]);
      y(l, k) = f(t[l], t[k]);
    }
  }
  // Loop over all evaluation points
  const std::vector<double>::size_type N = x.size();
  std::vector<double> res(N, 0.0);  // Result vector
  std::vector<double> wx(q);        // Weights $\cob{w_x^i}$
  std::vector<double> wy(q);        // Weights $\cob{w_y^i}$
  int nodex;  // Store index of tx in case of division by zero
  int nodey;  // Store index of ty in case of division by zero
  for (int k = 0; k < N; ++k) {
    double sx = 0;
    double sy = 0;
    bool nonodex = true;
    bool nonodey = true;

    for (int j = 0; j < q; ++j) {
      if (nonodex) {
        // In the case of division by zero
        if (std::abs(x[k][0] - t[j]) < 1e-9) {
          nonodex = false;
          nodex = j;       // Keep track of the index
          if (!nonodey) {  // Early termination in the case of second division
                           // by zero
            break;
          }
        } else {
          wx[j] = lambda[j] / (x[k][0] - t[j]);
          sx += wx[j];
        }
      }
      // Only enter if no division by zero occured
      if (nonodey) {
        // In the case of division by zero
        if (std::abs(x[k][0] - t[j]) < 1e-9) {
          nonodey = false;
          nodey = j;       // Keep track of the index
          if (!nonodex) {  // Early termination in the case of second division
                           // by zero
            break;
          }
        } else {
          wy[j] = lambda[j] / (x[k][1] - t[j]);
          sy += wy[j];
        }
      }
    }
    // The case of 2 division by zero
    if (!nonodex && !nonodey) {
      res[k] = y(nodex, nodey);
    }
    // The case of a division by zero along the x Axis
    else if (!nonodex) {
      for (int j = 0; j < q; ++j) {
        res[k] += y(nodex, j) * wy[j];
      }
      res[k] /= sy;
    }
    // The case of a division by zero along the y Axis
    else if (!nonodey) {
      for (int j = 0; j < q; ++j) {
        res[k] += y(j, nodey) * wx[j];
      }
      res[k] /= sx;
    }
    // No division by zero occurred
    else {
      for (int j = 0; j < q; ++j) {
        for (int i = 0; i < q; ++i) {
          res[k] += y(j, i) * wx[j] * wy[i];
        }
      }
      res[k] /= sy * sx;
    }
  }
  return res;
}
/* SAM_LISTING_END_2 */

/* SAM_LISTING_BEGIN_3 */
template <typename FUNCTOR>
std::vector<double> genChebInterpEval2D(unsigned int q, FUNCTOR f,
                                        Eigen::Vector2d a, Eigen::Vector2d b,
                                        const std::vector<Eigen::Vector2d> &x) {
  const std::vector<double>::size_type N = x.size();
  std::vector<double> res = std::vector<double>(N, 0.0);  // Result vector

  // Transformation from $\cob{\cintv{-1,1}^2}$ to $\cob{\cintv{a_1,b_1}\times\cintv{a_2,b_2}}$
  auto phif = [f, a, b](double x1, double x2) {
    const double tmp1 = 0.5 * ((b[0] - a[0]) * x1 + a[0] + b[0]);
    const double tmp2 = 0.5 * ((b[1] - a[1]) * x2 + a[1] + b[1]);
    return f(tmp1, tmp2);
  };
  // Inverse transformation $\cob{\Phibf^-1}$
  auto phiinv = [a, b](const Eigen::Vector2d &x) {
    const double x1 = 2. / (b[0] - a[0]) * (x[0] - 0.5 * (a[0] + b[0]));
    const double x2 = 2. / (b[1] - a[1]) * (x[1] - 0.5 * (a[1] + b[1]));
    return (Eigen::Vector2d() << x1, x2).finished();
  };
  // Vector of transformed evaluation points
  std::vector<Eigen::Vector2d> phiinvx(N);
  for (unsigned int k = 0; k < N; ++k) {
    phiinvx[k] = phiinv(x[k]); // Affine transformation
  }
  // Interpolation on $\cob{\cintv{-1,1}^2}$
  res = chebInterpEval2D(q, phif, phiinvx);
  return res;
}
/* SAM_LISTING_END_3 */

// Some functions for debugging   
template <typename FUNCTOR>
double errorestimate(unsigned int q, FUNCTOR f,
                     const std::vector<Eigen::Vector2d> &x,
                     const std::vector<double> &res) {
  const std::vector<double>::size_type N = x.size();
  double error = 0.0;
  for (int i = 0; i < N; ++i) {
    Eigen::Vector2d xi = x[i];
    error += std::pow(f(xi[0], xi[1]) - res[i], 2);
  }
  return error / N;
}
template <typename FUNCTOR>
double errorestimate(unsigned int q, FUNCTOR f, const std::vector<double> &x,
                     const std::vector<double> &res) {
  const std::vector<double>::size_type N = x.size();
  double error = 0.0;
  for (int i = 0; i < N; ++i) {
    error += std::pow(f(x[i]) - res[i], 2);
  }
  return error / N;
}
}  // namespace TensorProductChebIntp

#endif  // TENSORPRODCHEB_H_
