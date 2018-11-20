/**
 * \file abstract_parametrized_curve.hpp
 * \brief This file declares an abstract class for representing
 *        Parametrized curves.
 *
 *  This File is a part of the 2D-Parametric BEM package
 */

#ifndef PARAMETRIZEDCURVEHPP
#define PARAMETRIZEDCURVEHPP

#include <iterator>
#include <memory>
#include <utility>
#include <vector>
#include <iostream>

#include <Eigen/Dense>
#include "gradient_descent.hpp"
/**
 * \namespace parametricbem2d
 * \brief This namespace contains all the classes and functions for
 *        the 2D-Parametric BEM package.
 */
namespace parametricbem2d {
// Forward declaration of the class for using it in the typedef for PanelVector
class AbstractParametrizedCurve;

/**
 * This typedef is used for defining a vector of AbstractParametrizedCurve
 * pointers which can store any sequence of parametrized curves derived from the
 * abstract base class AbstractParametrizedCurve. The PanelVector is also used
 * to store the components of a mesh by using additional constraints. Smart
 * pointers are used to free up memory whenever the pointer is destroyed.
 */
using PanelVector = std::vector<std::shared_ptr<AbstractParametrizedCurve>>;

/**
 * \class AbstractParametrizedCurve
 * \brief This abstract class provides the interface for a parametric
 *        curve. The interface follows the one outlined in
 *        \f$\ref{cpp:crv}\f$ in the lecture material for
 *        Advanced Numerical Methods for CSE.
 */
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
  static std::pair<double, double> ParameterRange(void) {
    // Parameter range : [-1,1]
    return std::make_pair(-1., 1.);
  }

  /**
   * This function is used for accessing a point on the parametrized
   * curve with parameter "t" for the parametrization \f$\gamma\f$(t).
   * This is a pure virtual function which has to be implemented
   * in the inherited classes.
   *
   * @param t Double type argument for the parameter value.
   * @return A 2-D vector of type Eigen::Vector2d
   *         containing the coordinates of the
   *         parametrized point.
   */
  virtual Eigen::Vector2d operator()(double t) const = 0;

  /**
   * This function is used for retrieving the derivative \f$\dot{\gamma}\f$(t)
   * at a point on the parametrized curve \f$\gamma\f$(t).
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
   * This function is used for retrieving the double derivative
   * \f$\ddot{\gamma}\f$(t)
   * at a point on the parametrized curve \f$\gamma\f$(t).
   * This is a pure virtual function which has to be implemented
   * in the inherited classes.
   *
   * @param t Parameter value of Double type
   *          at which the double derivative for
   *          the parametrization is evaluated
   * @return A 2-D vector of type Eigen::Vector2d
   *         containing the parametrized double
   *         derivative at point 't'
   */
  virtual Eigen::Vector2d DoubleDerivative(double t) const = 0;

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
    double a, b;
    // Getting the parameter range
    std::tie(a, b) = ParameterRange();
    // Checking if the value is within parameter range
    return (t >= a && t <= b);
  }

  /**
   * This function is used for splitting a parametrized curve into several
   * self-similar part curves. It is useful to make a parametrized mesh as it
   * returns all the part-parametrizations stored in a sequence in a PanelVector
   * such that the end point of a part is the starting point of the next part.
   *
   * @param N Integer indicating the number of parts for split
   * @return PanelVector containing all the part parametrizations
   */
  virtual PanelVector split(unsigned int N) const = 0;

  /**
   * This function is used for splitting a parametrized curve into several
   * self-similar part curves. It is useful to make a parametrized mesh as it
   * returns all the part-parametrizations stored in a sequence in a PanelVector
   * such that the end point of a part is the starting point of the next part.
   *
   * @param N Integer indicating the number of parts for split
   * @return PanelVector containing all the part parametrizations
   */
  std::pair<Eigen::Vector2d, double>
  distanceTo(const AbstractParametrizedCurve &curve) const {
    // The vector to store the solution for minimum distance between two curves
    Eigen::Vector2d initial(2);
    initial << 0,0;
    Eigen::Vector2d interior_solution(2);
    // Finding an interior solution
    std::function<Eigen::Vector2d(Eigen::Vector2d)> gradient = [&] (Eigen::Vector2d x) {
      double s, t;
      s = x(0);
      t = x(1);
      return Eigen::Vector2d(
                (this->operator()(s) - curve(t)).dot(this->Derivative(s)),
                -(this->operator()(s) - curve(t)).dot(curve.Derivative(t)));
    };
    bool interior_soln_exists;
    std::tie(interior_solution,interior_soln_exists) = GradientDescent(gradient,initial);
    std::cout << "Interior solution " << interior_solution << std::endl;
    if (interior_soln_exists)
      std::cout << "interior solution exists " << std::endl;
    // Find minima on the boundary
    double guess = 0.;
    bool boundary;
    double solution;
    double min_dist;
    Eigen::Vector2d final_solution;
    // this curve(-1)
    /*std::function<double(double)> gradient1 = [&] (double t) {
      return -(this->operator()(-1) - curve(t)).dot(curve.Derivative(t));
    };
    std::tie(solution,boundary) = GradientDescent(gradient1,guess);
    if (boundary)
      solution = std::signbit(solution);
    double dist1 = (this->operator()(-1) - curve(solution)).norm();
    min_dist = dist1;
    final_solution = Eigen::Vector2d(-1,solution);

    // this curve(1)
    std::function<double(double)> gradient2 = [&] (double t) {
      return -(this->operator()(1) - curve(t)).dot(curve.Derivative(t));
    };
    std::tie(solution,boundary) = GradientDescent(gradient2,guess);
    if (boundary)
      solution = std::signbit(solution);
    double dist2 = (this->operator()(1) - curve(solution)).norm();
    if (dist2 < min_dist) {
      min_dist = dist2;
      final_solution = Eigen::Vector2d(1,solution);
    }

    // curve(-1)
    std::function<double(double)> gradient3 = [&] (double s) {
      return (this->operator()(s) - curve(-1)).dot(this->Derivative(s));
    };
    std::tie(solution,boundary) = GradientDescent(gradient3,guess);
    if (boundary)
      solution = std::signbit(solution);
    double dist3 = (this->operator()(solution) - curve(-1)).norm();
    if (dist3 < min_dist) {
      min_dist = dist3;
      final_solution = Eigen::Vector2d(solution,-1);
    }

    // curve(1)
    std::function<double(double)> gradient4 = [&] (double s) {
      return (this->operator()(s) - curve(1)).dot(this->Derivative(s));
    };
    std::tie(solution,boundary) = GradientDescent(gradient4,guess);
    if (boundary)
      solution = std::signbit(solution);
    double dist4 = (this->operator()(solution) - curve(1)).norm();
    if (dist4 < min_dist) {
      min_dist = dist4;
      final_solution = Eigen::Vector2d(solution,1);
    }

    */
    // gradient descent stopped on the boundary
    if (interior_soln_exists || true) {
      double s, t;
      s = interior_solution(0);
      t = interior_solution(1);
      double interior_distance = (this->operator()(s) - curve(t)).norm();
      if (interior_distance < min_dist || true) {
        min_dist = interior_distance;
        final_solution = interior_solution;
      }
    }

    return std::make_pair(final_solution,min_dist);
  }
}; // class AbstractParametrizedCurve
} // namespace parametricbem2d

#endif // PARAMETRIZEDCURVEHPP
