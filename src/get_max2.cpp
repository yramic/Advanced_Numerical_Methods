#include <math.h>

double get_max2(double xl,double xr,double yl,double yr){
  if(std::abs(xr-xl)>std::abs(yr-yl)) return std::abs(xr-xl);	
  else 	return std::abs(yr-yl);
};
