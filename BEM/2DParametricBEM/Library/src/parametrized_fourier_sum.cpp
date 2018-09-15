/**
 * \file parametrized_fourier_sum.cpp
 * \brief This file defines the class representing Fourier
 *        Sum based parametrization
 * @see parametrized_fourier_sum.hpp
 *
 *  This File is a part of the 2D-Parametric BEM package
 */

#include "abstract_parametrized_curve.hpp"
#include "parametrized_fourier_sum.hpp"

#include <math.h>
#include <assert.h>
#include <iostream>
#include <utility>
#include <vector>

#include <Eigen/Dense>



using CoefficientsList = typename ParametrizedFourierSum::CoefficientsList;

ParametrizedFourierSum::ParametrizedFourierSum(CoefficientsList cos_list,
                                               CoefficientsList sin_list)
                                              :cosine(cos_list),sine(sin_list)
                       {
                         // Checking consistency
                         assert(cosine.rows()==sine.rows() &&
                                cosine.cols()==sine.cols());
                       }

Eigen::Vector2d ParametrizedFourierSum::operator() (double t) const {
  assert(IsWithinParameterRange(t));
  int N = cosine.cols();
  Eigen::VectorXd cos_theta(N);
  Eigen::VectorXd sin_theta(N);
  for (int i = 0 ; i<N ; ++i) {
    // Computing cosine and sine values at the parameter value t
    cos_theta(i) = cos((i+1)*t);
    sin_theta(i) = sin((i+1)*t);
  }
  // Matrix multiplication to create the Fourier Sum (in vector form)
  Eigen::Vector2d point = cosine*cos_theta + sine*sin_theta;
  return point;
}

Eigen::Vector2d ParametrizedFourierSum::Derivative(double t) const {
  assert(IsWithinParameterRange(t));
  int N = cosine.cols();
  Eigen::VectorXd cos_theta_dot(N);
  Eigen::VectorXd sin_theta_dot(N);
  for (int i = 0 ; i<N ; ++i) {
    // Computing the derivatives for cosine and sine terms
    cos_theta_dot(i) = -(i+1)*sin((i+1)*t);
    sin_theta_dot(i) = (i+1)*cos((i+1)*t);
  }
  // Matrix multiplication to create the derivative of Fourier Sum (in vector form)
  Eigen::Vector2d derivative = cosine*cos_theta_dot + sine*sin_theta_dot;
  return derivative;
}
