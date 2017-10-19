#include "../include/kernel.hpp"
#include <cmath>
#include <limits>
#include <iostream>
#include "../include/cheby.hpp"

// kernel for 2D problem with grid
double Kernel2D::operator()(double x, double y)
{
    if(std::abs(x-y) > std::numeric_limits<double>::epsilon())
        return num_ / std::abs(x-y);
    else
        return 0.;
}
double Kernel2D::operator()(double x1, double y1, double x2, double y2){}
// kernel for 4D problem of 2 vectors
double Kernel4D::operator()(double x1, double y1, double x2, double y2)
{
    double lvl;
    lvl = std::sqrt(std::pow(x1 - x2, 2) + std::pow(y1 - y2, 2)); // should be ||x-y||
    return -(1./(2*M_PI))*(std::log(lvl));
}

// 4D Polynomial Kernel
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
double GlobalSmoothKernel::operator()(double x, double y)
{
    return std::cos(std::sqrt(std::abs(x-y)));
}

// Gauss Kernel
double GaussKernel::operator()(double x1, double y1, double x2, double y2)
{
    return std::exp(-(std::pow((std::sqrt(std::pow(x1 - x2, 2) + std::pow(y1 - y2, 2))),2)));
}
