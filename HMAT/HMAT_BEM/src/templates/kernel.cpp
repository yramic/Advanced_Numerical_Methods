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
#include "../../include/kernel.hpp"

#include <cmath>

// kernel for 2D problem of 2 vectors
double KernelGalerkin::operator()(double x1, double y1, double x2, double y2) {
  if (std::abs(x1 - x2) < std::numeric_limits<double>::epsilon() &&
      std::abs(y1 - y2) < std::numeric_limits<double>::epsilon())
    return 0;
  double lvl;
  lvl = std::pow(x1 - x2, 2) + std::pow(y1 - y2, 2);  // should be ||x-y||^2
  double lvl1 = -1. / (4 * M_PI) * std::log(lvl);
  return lvl1;
}

// kernel for 2D problem of 2 segments
double KernelGalerkin::operator()(const Eigen::Vector2d& a,
                                  const Eigen::Vector2d& b,
                                  const Eigen::Vector2d& c,
                                  const Eigen::Vector2d& d, double eta) {
  // TODO
}
