/**
 * \file abstract_parametrized_curve.hpp
 * \brief This file declares a pure abstract class for representing
 *        Parametrized curves. The interface follows the one outlined
 *        Section 1.4.3.4 in the lecture material for Advanced Numerical
 *        Methods for CSE.
 *
 *  This File is a part of the 2D-Parametric BEM package
 */

#ifndef PARAMETRIZEDCURVEHPP
#define PARAMETRIZEDCURVEHPP

#include <utility>
#include <Eigen/Dense>

class AbstractParametrizedCurve {
  public:
/**
 * This function is used for querying the parameter interval
 *
 * @return A std::pair<double,double> object containing
 *          the valid parameter range for parametrization
 */
  virtual std::pair<double,double> ParameterRange(void) const = 0;

/**
 * This function is used for accessing a point on the parametrized
 * curve with parameter "t" for the parametrization $\gamma$(t)
 *
 * @param t Double type argument for the parameter value/
 * @return A 2-D vector of type Eigen::Vector2d
 *         containing the coordinates of the
 *         parametrized point.
 */
  virtual Eigen::Vector2d operator() (double t) const = 0;

/**
 * This function is used for retrieving the derivative \dot{$\gamma$(t)}
 * at a point on the parametrized curve
 *
 * @param t Parameter value of Double type
 *          at which the derivative for
 *          parametrization is evaluated
 * @return A 2-D vector of type Eigen::Vector2d
 *         containing the parametrized
 *         derivative at point 't'
 */
  virtual Eigen::Vector2d Derivative(double t) const = 0;
};

#endif //PARAMETRIZEDCURVEHPP
