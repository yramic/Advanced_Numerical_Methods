/**
 * \file parametrized_fourier_sum.cpp
 * \brief This file defines the class representing Fourier
 *        Sum based parametrization
 * @see parametrized_fourier_sum.hpp
 *
 *  This File is a part of the 2D-Parametric BEM package
 */

#include "parametrized_fourier_sum.hpp"

#include <math.h>
#include <assert.h>
#include <iostream>
#include <utility>
#include <vector>

#include <Eigen/Dense>


namespace parametricbem2d {
  using CoefficientsList = typename ParametrizedFourierSum::CoefficientsList;

  ParametrizedFourierSum::ParametrizedFourierSum(CoefficientsList cos_list,
                                                 CoefficientsList sin_list,
                                                 double tmin,
                                                 double tmax)
                                                :cosine_(cos_list),
                                                 sine_(sin_list),
                                                 tmin_(tmin),
                                                 tmax_(tmax)
                         {
                           // Checking consistency
                           assert(cosine_.cols()==sine_.cols());
                         }

  Eigen::Vector2d ParametrizedFourierSum::operator() (double t) const {
    assert(IsWithinParameterRange(t));
    t = t*(tmax_-tmin_)/2+(tmax_+tmin_)/2;
    int N = cosine_.cols();
    Eigen::VectorXd cos_theta(N);
    Eigen::VectorXd sin_theta(N);
    for (int i = 0 ; i<N ; ++i) {
      // Computing cosine and sine values at the parameter value t
      cos_theta(i) = cos((i+1)*t);
      sin_theta(i) = sin((i+1)*t);
    }
    // Matrix multiplication to create the Fourier Sum (in vector form)
    Eigen::Vector2d point = cosine_*cos_theta + sine_*sin_theta;
    return point;
  }

  Eigen::Vector2d ParametrizedFourierSum::Derivative(double t) const {
    assert(IsWithinParameterRange(t));
    t = t*(tmax_-tmin_)/2+(tmax_+tmin_)/2;
    int N = cosine_.cols();
    Eigen::VectorXd cos_theta_dot(N);
    Eigen::VectorXd sin_theta_dot(N);
    for (int i = 0 ; i<N ; ++i) {
      // Computing the derivatives for cosine and sine terms
      cos_theta_dot(i) = -(i+1)*sin((i+1)*t);
      sin_theta_dot(i) = (i+1)*cos((i+1)*t);
    }
    // Matrix multiplication to create the derivative of Fourier Sum (in vector form)
    Eigen::Vector2d derivative = cosine_*cos_theta_dot + sine_*sin_theta_dot;
    return derivative;
  }

  PanelVector ParametrizedFourierSum::split(unsigned int N) const {
    PanelVector parametrization_parts;
    for (int i = 0 ; i < N ; ++i) {
      double tmin = tmin_+ i*(tmax_-tmin_)/N;
      double tmax = tmin_+ (i+1)*(tmax_-tmin_)/N;
      ParametrizedFourierSum part(cosine_,sine_,tmin,tmax);
      parametrization_parts.push_back(&part);
    }
    return parametrization_parts;
  }
} // namespace parametricbem2d
