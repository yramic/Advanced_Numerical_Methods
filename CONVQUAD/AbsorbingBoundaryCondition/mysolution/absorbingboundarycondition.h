
#ifndef ABC_H_
#define ABC_H_


#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <Eigen/SparseLU>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <unsupported/Eigen/FFT>
using namespace Eigen;
using namespace std;

/* @brief Compute the convolution quadrature weights for Laplace transform F
 * \param F     Template function for the Laplace transform
 * \param delta Template function determining the multistep method
 * \param tau   Time step size  
 * \param M     Number of time steps
 * \\return convolution quadrature weights
 */
/* SAM_LISTING_BEGIN_0 */
template <typename FFUNC, typename DFUNC>
VectorXd cqweights_by_dft(const FFUNC& F, const DFUNC& delta, double tau,
                          size_t M) {
  Eigen::VectorXcd w = Eigen::VectorXd::Zero(M + 1);
