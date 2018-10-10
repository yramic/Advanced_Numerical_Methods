/**
 * \file parametrized_line.cpp
 * \brief This file defines the class representing parametrization of
 *        a line segment in 2D.
 * @see parametrized_line.hpp
 *
 *  This File is a part of the 2D-Parametric BEM package
 */

#include "parametrized_line.hpp"

#include <math.h>
#include <assert.h>
#include <iostream>
#include <utility>

#include <Eigen/Dense>

namespace parametricbem2d {
  using Point = typename ParametrizedLine::Point;

  ParametrizedLine::ParametrizedLine(Point first, Point second) : start_(first),end_(second) {}

  Point ParametrizedLine::operator() (double t) const {
    assert(IsWithinParameterRange(t));
    double x1(start_(0)),y1(start_(1)),x2(end_(0)),y2(end_(1));
    // Linear interpolation of x & y coordinates based on parameter t
    Eigen::Vector2d point(t*(x2-x1)/2 + (x2+x1)/2
                         ,t*(y2-y1)/2 + (y2+y1)/2);
    return point;
  }

  Eigen::Vector2d ParametrizedLine::Derivative(double t) const {
    assert(IsWithinParameterRange(t));
    double x1(start_(0)),y1(start_(1)),x2(end_(0)),y2(end_(1));
    // Derivative of the linear interpolation used in the function operator()
    Eigen::Vector2d derivative((x2-x1)/2,(y2-y1)/2);
    return derivative;
  }

  Eigen::Vector2d ParametrizedLine::DoubleDerivative(double t) const {
    assert(IsWithinParameterRange(t));
    // Double Derivative of the linear interpolation used in function operator()
    Eigen::Vector2d double_derivative(0,0);
    return double_derivative;
  }

  PanelVector ParametrizedLine::split(unsigned int N) const {
    PanelVector parametrization_parts;
    for (int i = 0 ; i < N ; ++i) {
      double tmin,tmax;
      std::tie(tmin,tmax) = ParameterRange();
      double t1 = tmin + i*(tmax-tmin)/N;
      double t2 = tmin + (i+1)*(tmax-tmin)/N;
      Point first = this->operator()(t1);
      Point second = this->operator()(t2);
      parametrization_parts.push_back(std::make_shared<ParametrizedLine>
                                      (first,second));
    }
    return parametrization_parts;
  }
} // namespace parametricbem2d
