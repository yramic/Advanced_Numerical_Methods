/***********************************************************************
 *                                                                     *
 * Code for Course "Advanced Numerical Methods for CSE"                *
 * (Prof. Dr. R. Hiptmair)                                             *
 * Author: Daniele Casati                                              *
 * Date: 11/2017                                                       *
 * (C) Seminar for Applied Mathematics, ETH Zurich                     *
 * This code can be freely used for non-commercial purposes as long    *
 * as this header is left intact.                                      *
 ***********************************************************************/
#include <cmath>

// get max edge of the cluster
double get_max(double xl, double xr, double yl, double yr) {
  if (std::abs(xr - xl) > std::abs(yr - yl))
    return std::abs(xr - xl);
  else
    return std::abs(yr - yl);
}

// get distance from diagonal
double get_min(double xl, double xr, double yl, double yr) {
  double dist = std::abs(xl - yl);
  double dist2 = std::abs(xl - yr);
  double dist3 = std::abs(xr - yl);
  double dist4 = std::abs(xr - yr);

  if (dist2 < dist) dist = dist2;
  if (dist3 < dist) dist = dist3;
  if (dist4 < dist) dist = dist4;

  return dist;
}

// check the admissibility of a bounding box
bool is_admissible(double xl, double xr, double yl, double yr, double eta) {
  return get_max(xl, xr, yl, yr) <= eta * get_min(xl, xr, yl, yr);
}
