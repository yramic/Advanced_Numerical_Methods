/*
  This hpp file declares the pure abstract class for representing parametrized
  Curves according to the interface outlined in Section 1.4.3.4.
*/

#ifndef PARAMETRIZEDCURVEHPP
#define PARAMETRIZEDCURVEHPP

#include <utility>
#include <Eigen/Dense>

class AbstractParametrizedCurve {
  public:
  /* This function is used for querying the parameter interval
  *
  * Parameters : None
  *
  * Return : std::pair<double,double>
  *          The valid parameter range for parametrization
  */
  virtual std::pair<double,double> ParameterRange(void) const = 0;

  /* This function is used for accessing a point on the parametrized
  * curve with parameter "t" $\gamma$(t)
  *
  * Parameters :
  *
  * t : double
  *     Parameter value at which the parametrization point
  *     is evaluated
  *
  * Return : Eigen::Vector2d
  *          2-D vector containing the coordinates of the
  *          parametrized point
  */
  virtual Eigen::Vector2d operator() (double t) const = 0;

  /* This function is used for retrieving the derivative \dot{$\gamma$(t)}
  * at a point on the parametrized curve
  *
  * Parameters :
  *
  * t : double
  *     Parameter value at which the derivative for
  *     parametrization is evaluated
  *
  * Return : Eigen::Vector2d
  *          2-D vector containing the parametrized
  *          derivative at point 't'
  */
  virtual Eigen::Vector2d Derivative(double t) const = 0;
};

#endif //PARAMETRIZEDCURVEHPP
