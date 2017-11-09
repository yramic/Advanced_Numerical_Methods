#include "../../include/kernel.hpp"
#include <cmath>
#include <limits>

// Kernel functor \f$\log{\left|x-y\right|}\f$ if $x != y$, else 0
double KernelLog::operator()(double x, double y)
{
    if(std::abs(x-y) > std::numeric_limits<double>::epsilon())
        return num_ * std::log(std::abs(x-y));
    else
        return 0.;
}

// Kernel functor $\frac{num}{|x-y|}$ if $x != y$, else 0
double KernelInvDistance::operator()(double x, double y)
{
    // TODO
    return 0.;
}

// Kernel functor \f$\frac{\cos(C\left|x-y\right|)}{\left|x-y\right|}\f$ if $x != y$, else 0
double KernelCosine::operator()(double x, double y)
{
    // TODO
    return 0.;
}
