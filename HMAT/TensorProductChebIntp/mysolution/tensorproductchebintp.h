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
    t[k] =
        std::cos((2.0 * (k + 1) - 1.0) / (2 * q) * M_PI);  // \prbeqref{eq:chn}
    lambda[k] =
        std::pow(2, q - 1) / q * sgn *
        std::sin((2.0 * (k + 1) - 1.0) / (2 * q) * M_PI);  // \prbeqref{eq:bwf}
    y[k] = f(t[k]);                                        // $\cob{y_k}$
  }
  // Loop over all evaluation points $\cob{x_i}$
  const std::vector<double>::size_type N = x.size();
  std::vector<double> res(N, 0.0);  // Result vector
  for (int k = 0; k < N; ++k) {
    // **********************************************************************
    // Your Solution here -> TASK 2-3.c!
    // First we need to differentiate between two cases: 
    // Case one: x[i] = t[i] in this case the denominator would turn 0 which is a problem!
    // Case two: x[i] != t[i] in this case we can actually do the interpolation!

    bool nonode {true}; // Check if x[i] != t[i]
    double denom {0}; // Initialization of the denominator
    for (unsigned int j {0}; j < q; j++) {
      // Condition 1:
      if (x[k] == t[j]) {
        res[k] = y[j];
        nonode = false;
        break; // Jump out of the loop!
      }
      // Now we need to implement Eq. (2.3.4.16) from the lecture slides!
      const double tmp = lambda[j] /(x[k] - t[j]);
      res[k] += tmp * y[j];
      denom += tmp;
    }
    // Now we are out of the j loop and we need to devide everything with the denominator
    // This is because the formula looks something like this: denom * (denom)^-1 ! 
    // But this should be done only if x[i] != t[i] otherwise we sould break again!
    if (nonode) {
      res[k] /= denom;
    }
    // **********************************************************************
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
    t[k] = std::cos((2.0 * (k + 1) - 1.0) / (2 * q) * M_PI);
    lambda[k] = std::pow(2, q - 1) / q * sgn *
                std::sin((2.0 * (k + 1) - 1.0) / (2 * q) * M_PI);
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

    // **********************************************************************
    // Your Solution here for Problem (2-3.d)
    // Here we have the same problem but now in 2D, hence we need to check,
    // if x != t and y != t as well!
    for (int j{0}; j < q; ++j) {
      if (x[k][0] == t[j]) {
        nonodex = false;
        nodex = j;
        if (!nonodey) {
          break;
        }
      } else {
        wx[j] = lambda[j] / (x[k][0] - t[j]);
        sx += wx[j];
      }

      if (x[k][1] == t[j]) {
        nonodey = false;
        nodey = j;
        if (!nonodex) {
          break;
        }
      } else {
        wy[j] = lambda[j] / (x[k][1] -t[j]);
        sy += wy[j];
      }
    }

    // Now we need to differentiate between 3 cases
    // Case 1: y = t[j] && x = t[j]
    if (!nonodex && !nonodey) {
      res[k] = y(nodex, nodey);
    } 
    // Case 2: y = t[j] && x != t[j]
    else if (nonodex && !nonodey) {
      for (int j{0}; j < q; ++j) {
        res[k] = wx[j] * y(j, nodey);
      }
      res[k] /= sx;
    }
    // Case 3: y != t[j] && x = t[j]
    else if (nonodex && !nonodey) {
      for (int j{0}; j < q; ++j) {
        res[k] = wy[j] * y(nodex, j);
      }
      res[k] /= sy;
    }
    // Case 4: y!= t[j] && x != t[j]
    else {
      for (int i{0}; i < q; ++i) {
        for (int j{0}; j < q; ++j) {
          res[k] = wx[i] * wy[j] * y(i,j);
        }
      }
      res[k] /= (sx *sy);
    }

    // **********************************************************************
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

  // **********************************************************************
  // Your Solution here for Problem 2-3e:
  // In this task I want to make an affine transformation from [-1,1] to [a,b]
  // But the problem is still in 2D Hence I have to work with matrices!
 
  // Define the transformation matrix
  Eigen::MatrixXd m = Eigen::MatrixXd::Identity(2, 2);
  m(0, 0) *= (2.0 / (b.x() - a.x()));
  m(1, 1) *= (2.0 / (b.y() - a.y()));

  // Define the translation vector
  Eigen::Vector2d v(0.5 * (b.x() + a.x()), 0.5 * (b.y() + a.y()));

  // Transform and evaluate for each point
  for (int i = 0; i < N; ++i) {
    // Apply the affine transformation to the input point 'x[i]'
    Eigen::Vector2d x_transformed = m * (x[i] - v);

    // Evaluate the function 'f' on the transformed point
    double result = f(x_transformed[0], x_transformed[1]);

    // Store the result
    res[i] = result;
  }

  // **********************************************************************
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
