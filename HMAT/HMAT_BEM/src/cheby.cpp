/***********************************************************************
 *                                                                     *
 * Code for Course "Advanced Numerical Methods for CSE"                *
 * (Prof. Dr. R. Hiptmair)                                             *
 * Author: Daniele Casati                                              *
 * Date: 11/2017                                                       *
 * (C) Seminar for Applied Mathematics, ETH Zurich                     *
 * This code can be freely used for non-commercial purposes as long    *
 * as this header is left intact.                                      *
 ***********************************************************************/
#include "../include/cheby.hpp"

#include <Eigen/Dense>
#include <cmath>

// constructor
Cheby::Cheby(double xl, double xr, unsigned deg)
    : xl_(xl),
      xr_(xr),
      deg_(deg),
      tk_(Eigen::VectorXd::Zero(deg + 1)),
      wk_(Eigen::VectorXd::Zero(deg + 1)) {
  setNodes();
  setWghts();
}

// compute Chebyshew nodes on domain [xl,xr]
void Cheby::setNodes() {
  tk_ = Eigen::VectorXd::Zero(deg_ + 1);
  for (unsigned j = 0; j <= deg_; ++j) {
    tk_(j) = (xl_ + xr_) / 2 +
             (xr_ - xl_) / 2 * cos((2. * j + 1) / (2. * (deg_ + 1.)) * M_PI);
  }
}

// compute weights of Lagrange polynomial
void Cheby::setWghts() {
  wk_ = Eigen::VectorXd::Zero(deg_ + 1);
  for (unsigned j = 0; j <= deg_; ++j) {
    double hc = 1.;
    for (unsigned k = 0; k < j; ++k) hc *= tk_[j] - tk_[k];
    // Skip ``k == j''
    for (unsigned k = j + 1; k <= deg_; ++k) hc *= tk_[j] - tk_[k];

    wk_(j) = 1. / hc;
  }
}
