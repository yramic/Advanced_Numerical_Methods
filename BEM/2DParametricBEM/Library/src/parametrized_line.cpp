/*
 * This cpp file defines the class representing
 * parametrization of a line in 2-D
 * (See parametrized_line.hpp)
 */

#include "abstract_parametrized_curve.hpp"
#include "parametrized_line.hpp"

#include <math.h>
#include <assert.h>
#include <iostream>
#include <utility>

#include <Eigen/Dense>


using Point = typename ParametrizedLine::Point;

ParametrizedLine::ParametrizedLine(Point first, Point second) : start(first),end(second) {}

std::pair<double,double> ParametrizedLine::ParameterRange(void) const {
  // Parameter range : [-1,1]
  return std::make_pair(-1.,1.);
}

Point ParametrizedLine::operator() (double t) const {
  assert(t>=-1 && t<=1);
  double x1(start(0)),y1(start(1)),x2(end(0)),y2(end(1));
  // Linear interpolation of x & y coordinates based on parameter t
  Eigen::Vector2d point(t*(x2-x1)/2 + (x2+x1)/2
                       ,t*(y2-y1)/2 + (y2+y1)/2);
  return point;
}

Eigen::Vector2d ParametrizedLine::Derivative(double t) const {
  assert(t>=-1 && t<=1);
  double x1(start(0)),y1(start(1)),x2(end(0)),y2(end(1));
  // Derivative of the linear interpolation used in the function operator()
  Eigen::Vector2d derivative((x2-x1)/2,(y2-y1)/2);
  return derivative;
}
