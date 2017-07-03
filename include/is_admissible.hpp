#ifndef IS_ADMISSIBLE_HPP
#define IS_ADMISSIBLE_HPP


// primitive functions used to check whether a cluster is admissible

double get_max(double xl, double xr, double yl, double yr);
double get_min(double xl, double xr, double yl, double yr);
bool is_admissible(double xl, double xr, double yl, double yr, double eta);

#endif // IS_ADMISSIBLE_HPP
