#include <math.h>
#include "../include/comparison.hpp"

bool is_admissible(double xl,double xr,double yl,double yr, double eta){
  return get_max2(xl,xr,yl,yr)<=(eta*get_min(xl,xr,yl,yr));
};
