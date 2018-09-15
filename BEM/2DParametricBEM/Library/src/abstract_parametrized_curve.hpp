/**
 * \file abstract_parametrized_curve.hpp
 * \brief This file declares an abstract class for representing
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
 * This function is used for querying the parameter interval.
 * The standard parameter interval [-1,1] is used and it can't
 * be overriden in the inherited classes as the function is non
 * virtual. The function is declared static as it is independent
 * of the concrete object.
 *
 * @return A std::pair<double,double> object containing
 *          the valid parameter range for parametrization
 *          that is [-1,1]
 */
 static std::pair<double,double> ParameterRange(void) {
   // Parameter range : [-1,1]
   return std::make_pair(-1.,1.);
 }

/**
 * This function is used for accessing a point on the parametrized
 * curve with parameter "t" for the parametrization $\gamma$(t).
 * This is a pure virtual function which has to be implemented
 * in the inherited classes.
 *
 * @param t Double type argument for the parameter value/
 * @return A 2-D vector of type Eigen::Vector2d
 *         containing the coordinates of the
 *         parametrized point.
 */
  virtual Eigen::Vector2d operator() (double t) const = 0;

/**
 * This function is used for retrieving the derivative \dot{$\gamma$(t)}
 * at a point on the parametrized curve.
 * This is a pure virtual function which has to be implemented
 * in the inherited classes.
 *
 * @param t Parameter value of Double type
 *          at which the derivative for
 *          parametrization is evaluated
 * @return A 2-D vector of type Eigen::Vector2d
 *         containing the parametrized
 *         derivative at point 't'
 */
  virtual Eigen::Vector2d Derivative(double t) const = 0;

/**
 * This function is used for checking whether a value t is within the
 * valid parameter range. This function is non virtual to prevent it
 * from being overriden as the parameter interval is fixed. It is
 * declared static because the check is independent of the concrete
 * implementation.
 *
 * @param t The value to be checked
 * @return boolean indicating result of the performed check
 */
   static bool IsWithinParameterRange(double t) {
     double a,b;
     std::tie(a,b) = ParameterRange();
     return ( t>=a && t<=b );
   }
};

#endif //PARAMETRIZEDCURVEHPP
