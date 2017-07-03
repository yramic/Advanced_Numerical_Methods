#include "../include/kernel.hpp"
#include <cmath>
#include <limits>

Kernel::Kernel( double num ):
    num_( num )
{ }

double Kernel::operator()( double x, double y ) {

    if(std::abs(x-y) > std::numeric_limits<double>::epsilon())
        return num_/std::abs(x-y);
    else
        return 0.;
}
