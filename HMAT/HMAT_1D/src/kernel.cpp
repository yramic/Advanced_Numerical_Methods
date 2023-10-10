#include "../include/kernel.hpp"

#include <cmath>
#include <limits>

// Kernel functor \f$\log{\left|x-y\right|}\f$ if $x != y$, else 0
double KernelLog::operator()(double x, double y) {
  if (std::abs(x - y) > std::numeric_limits<double>::epsilon())
    return num_ * std::log(std::abs(x - y));
  else
    return 0.;
}
