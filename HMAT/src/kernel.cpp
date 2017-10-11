#include "../include/kernel.hpp"
#include <cmath>
#include <limits>


double Kernel::operator()(double x, double y)
{
    if(std::abs(x-y) > std::numeric_limits<double>::epsilon())
        return num_ / std::abs(x-y);
    else
        return 0.;
}

double Kernel::operator()(double x1, double y1, double x2, double y2)
{
    double lvl;
    lvl = std::sqrt(std::pow(x1 - x2, 2) + std::pow(y1 - y2, 2)); // should be ||x-y||
    return -(1./(2*M_PI))*(std::log(lvl));
}
