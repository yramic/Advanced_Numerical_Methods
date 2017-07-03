#include "../include/is_admissible.hpp"
#include <cmath>


double get_max(double xl, double xr, double yl, double yr)
{
    if(std::abs(xr-xl)>std::abs(yr-yl))
        return std::abs(xr-xl);
    else
        return std::abs(yr-yl);
}


double get_min(double xl, double xr, double yl, double yr)
{
    double min4 = std::abs(xl-yl);
    double k2 = std::abs(xl-yr);
    double k3 = std::abs(xr-yl);
    double k4 = std::abs(xr-yr);

    if(k2 < min4)
        min4 = k2;
    if(k3 < min4)
        min4 = k3;
    if(k4 < min4)
        min4 = k4;

    return min4;
}


bool is_admissible(double xl, double xr, double yl, double yr, double eta)
{
    return get_max(xl,xr,yl,yr) <= eta*get_min(xl,xr,yl,yr);
}
