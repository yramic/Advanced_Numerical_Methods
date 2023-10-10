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
#include <cmath>
#include <limits>

// kernel for 2D problem of 2 vectors
double KernelGalerkin::operator()(double x1, double y1, double x2, double y2)
{
    if(std::abs(x1-x2) < std::numeric_limits<double>::epsilon() &&
       std::abs(y1-y2) < std::numeric_limits<double>::epsilon()) return 0;
    double lvl;
    lvl = std::pow(x1 - x2, 2) + std::pow(y1 - y2, 2); // should be ||x-y||^2
    double lvl1 = -1./(4*M_PI)*std::log(lvl);
    return num_ * lvl1;
}

// 2D Polynomial Kernel
double KernelPolynomial::operator()(double x1, double y1, double x2, double y2)
{
    return num_ * x1*x2*y1*y2;
}

// Constant Kernel
double KernelConstant::operator()(double x1, double y1, double x2, double y2)
{
    return num_;
}

// Global Smooth Kernel
double KernelGlobalSmooth::operator()(double x1, double y1, double x2, double y2)
{
    return num_ * std::cos(std::sqrt(std::pow(x1 - x2, 2) + std::pow(y1 - y2, 2)));
}

// Gauss Kernel
double KernelGauss::operator()(double x1, double y1, double x2, double y2)
{
    int t1 = std::pow(std::sqrt(std::pow(x1 - x2, 2) + std::pow(y1 - y2, 2)),2);
    int t  = std::exp(-(t1));
    return num_ * std::exp(-(std::pow((std::sqrt(std::pow(x1 - x2, 2) + std::pow(y1 - y2, 2))),2)));
}
