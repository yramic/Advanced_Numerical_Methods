/**
 * \file parametrized_circular_arc.cpp
 * \brief This file defines the class for representing parametrization
 *        of a circular arc.
 *
 *  This File is a part of the 2D-Parametric BEM package
 */

#include "parametrized_circular_arc.hpp"

#include <utility>

#include <Eigen/Dense>

namespace parametricbem2d {

  ParametrizedCircularArc::ParametrizedCircularArc(Eigen::Vector2d center,
                                                   double r,
                                                   double phi_start,
                                                   double phi_end):
                                                   center_(center),
                                                   radius_(r),
                                                   phi_start_(phi_start),
                                                   phi_end_(phi_end) {}

  Eigen::Vector2d ParametrizedCircularArc::operator() (double t) const {
    assert(IsWithinParameterRange(t));
    // Parametrization using polar coordinates based on parameter t
    double mean = (phi_start_ + phi_end_)/2.;
    double difference = (phi_end_ - phi_start_)/2.;
    Eigen::Vector2d point(center_(0) + radius_*cos(t*difference + mean),
                          center_(1) + radius_*sin(t*difference + mean));
    return point;
  }

  Eigen::Vector2d ParametrizedCircularArc::Derivative(double t) const {
    assert(IsWithinParameterRange(t));
    // Derivative of the polar coordinaties used in the function operator()
    double mean = (phi_start_ + phi_end_)/2.;
    double difference = (phi_end_ - phi_start_)/2.;
    Eigen::Vector2d derivative(-radius_*difference*sin(t*difference + mean),
                                radius_*difference*cos(t*difference + mean));
    return derivative;
  }

  PanelVector ParametrizedCircularArc::split(unsigned int N) const {
    PanelVector parametrization_parts;
    for (int i = 0 ; i < N ; ++i) {
      double phi_start = phi_start_ + i*(phi_end_-phi_start_)/N;
      double phi_end = phi_start_ + (i+1)*(phi_end_-phi_start_)/N;
      ParametrizedCircularArc part(center_,radius_,phi_start,phi_end);
      parametrization_parts.push_back(&part);
    }
    return parametrization_parts;
  }

} // namespace parametricbem2d
