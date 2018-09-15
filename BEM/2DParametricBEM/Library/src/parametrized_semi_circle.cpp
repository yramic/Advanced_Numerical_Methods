/**
 * \file parametrized_semi_circle.cpp
 * \brief This file defines the class representing parametrization
 *        of a semi circle.
 * @see parametrized_semi_circle.hpp
 *
 *  This File is a part of the 2D-Parametric BEM package
 */

#include "abstract_parametrized_curve.hpp"
#include "parametrized_semi_circle.hpp"

#include <math.h>
#include <assert.h>
#include <iostream>
#include <utility>

#include <Eigen/Dense>


#define _USE_MATH_DEFINES //for Pi

ParametrizedSemiCircle::ParametrizedSemiCircle(double r) : radius(r) {}

std::pair<double,double> ParametrizedSemiCircle::ParameterRange(void) const {
  // Parameter range : [-1,1]
  return std::make_pair(-1.,1.);
}

Eigen::Vector2d ParametrizedSemiCircle::operator() (double t) const {
  assert(t>=-1 && t<=1);
  // Parametrization using polar coordinates based on parameter t
  Eigen::Vector2d point(radius*cos(M_PI*t/2.),radius*sin(M_PI*t/2.));
  return point;
}

Eigen::Vector2d ParametrizedSemiCircle::Derivative(double t) const {
  assert(t>=-1 && t<=1);
  // Derivative of the polar coordinaties used in the function operator()
  Eigen::Vector2d derivative(-radius*M_PI*sin(M_PI*t/2.)/2.,M_PI*radius*cos(M_PI*t/2.)/2.);
  return derivative;
}
