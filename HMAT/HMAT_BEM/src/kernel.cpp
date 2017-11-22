/***********************************************************************
 *                                                                     *
 * Code for Course "Advanced Numerical Methods for CSE"                *
 * (Prof. Dr. R. Hiptmair)                                             *
 * Author: Ioannis Magkanaris                                          *
 * Date: 11/2017                                                       *
 * (C) Seminar for Applied Mathematics, ETH Zurich                     *
 * This code can be freely used for non-commercial purposes as long    *
 * as this header is left intact.                                      *
 ***********************************************************************/
#include "../include/kernel.hpp"
#include "../BEM/CppHilbert/Library/source/buildV.hpp"

// kernel for 2D problem of 2 vectors
double KernelGalerkin::operator()(const Eigen::Vector2d& a, const Eigen::Vector2d& b,
                                  const Eigen::Vector2d& c, const Eigen::Vector2d& d, double eta)
{
    return computeVij(a, b, c, d);
}
