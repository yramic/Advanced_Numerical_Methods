#include "../include/kernel.hpp"
#include <cmath>
#include <limits>
#include <iostream>
#include "../include/cheby.hpp"


// kernel for 2D problem of 2 vectors
double KernelGalerkin::operator()(double x1, double y1, double x2, double y2)
{
    if (x1== x2 || y1==y2) return 0;
    double lvl;
    lvl = std::sqrt(std::pow(x1 - x2, 2) + std::pow(y1 - y2, 2)); // should be ||x-y||
    double lvl1 = -(1./(2*M_PI))*(std::log(lvl));
    return lvl1;
}

// 2D Singular Kernel
double SingularKernel::operator()(double x1, double y1, double x2, double y2)
{
    if (x1== x2 || y1==y2) return 0;
    double lvl;
    lvl = std::sqrt(std::pow(x1 - x2, 2) + std::pow(y1 - y2, 2)); // should be ||x-y||
    double lvl1 = std::log(lvl);
    return lvl1;
}

// 2D Singular Kernel
double SingularKernelf::operator()(double x1, double y1, double x2, double y2)
{
    if (x1== x2 || y1==y2) return 0;
    double lvl;
    lvl = std::sqrt(std::pow(x1 - x2, 2) + std::pow(y1 - y2, 2)); // should be ||x-y||
    double lvl1 = std::log(1/lvl);
    return lvl1;
}
// 2D Polynomial Kernel
double PolynomialKernel::operator()(double x1, double y1, double x2, double y2)
{
    return x1*x2*y1*y2;
}

//Constant Kernel
double ConstantKernel::operator()(double x1, double y1, double x2, double y2)
{
    return num_;
}

// Global Smooth Kernel
double GlobalSmoothKernel::operator()(double x1, double y1, double x2, double y2)
{
    return std::cos(std::sqrt(std::pow(x1 - x2, 2) + std::pow(y1 - y2, 2)));
}


// Gauss Kernel
double GaussKernel::operator()(double x1, double y1, double x2, double y2)
{
    int t1 = std::pow(std::sqrt(std::pow(x1 - x2, 2) + std::pow(y1 - y2, 2)),2);
    int t = std::exp(-(t1));
    return std::exp(-(std::pow((std::sqrt(std::pow(x1 - x2, 2) + std::pow(y1 - y2, 2))),2)));
}
