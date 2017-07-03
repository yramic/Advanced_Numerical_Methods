#include <math.h>

double get_min(double xl, double xr, double yl, double yr){
  double min4=std::abs(xl-yl);
  double k2=std::abs(xl-yr);
  double k3=std::abs(xr-yl);
  double k4=std::abs(xr-yr);
  if(k2<min4)
    min4=k2;
  if(k3<min4)
    min4=k3;
  if(k4<min4)
    min4=k4;
  return min4;
};
