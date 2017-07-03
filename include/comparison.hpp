// primitive functions which are used to check whether a cluster is admissible
#ifndef COMPARISON_H
#define COMPARISON_H

#include <math.h>

double get_max2(double xl,double xr,double yl,double yr);
double get_min(double xl, double xr, double yl, double yr);
bool is_admissible(double xl,double xr,double yl,double yr, double eta);

#endif 
