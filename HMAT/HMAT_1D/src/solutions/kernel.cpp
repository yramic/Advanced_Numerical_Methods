#include "../../include/kernel.hpp"
#include <cmath>
#include <limits>

// Kernel functor $\frac{num}{|x-y|}$ if $x != y$, else 0
/* SAM_LISTING_BEGIN_0 */
double KernelInvDistance::operator()(double x, double y)
{
    if(std::abs(x-y) > std::numeric_limits<double>::epsilon())
        return num_ / std::abs(x-y);
    else
        return 0.;
}
/* SAM_LISTING_END_0 */

// Kernel functor \f$\frac{\cos(C\left|x-y\right|)}{\left|x-y\right|}\f$ if $x != y$, else 0
/* SAM_LISTING_BEGIN_1 */
double KernelCosine::operator()(double x, double y)
{
    if(std::abs(x-y) > std::numeric_limits<double>::epsilon())
        return std::cos(num_*std::abs(x-y)) / std::abs(x-y);
    else
        return 0.;
}
/* SAM_LISTING_END_1 */
