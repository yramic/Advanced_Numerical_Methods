#include "../../include/kernel.hpp"
#include <cmath>
#include <limits>

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
