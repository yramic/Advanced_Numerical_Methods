/**
 * \file parametrized_semi_circle.hpp
 * \brief This file declares the class for representing parametrization
 *        of a semi circle. It is centered at (0,0) and lies in
 *        positive x half plane. This class inherits from the
 *        Abstract base class representing parametrized curves
 * @see abstract_parametrized_curve.hpp
 *
 *  This File is a part of the 2D-Parametric BEM package
 */

#ifndef PARAMETRIZEDSEMICIRCLEHPP
#define PARAMETRIZEDSEMICIRCLEHPP

#include <utility>

#include <Eigen/Dense>

#include "abstract_parametrized_curve.hpp"

class ParametrizedSemiCircle : public AbstractParametrizedCurve {
  public:
    /**
     * Constructor with specified radius; default value = 1.
     */
    ParametrizedSemiCircle(double = 1.);

    /**
     * @see AbstractParametrizedCurve::ParameterRange()
     */
    std::pair<double,double> ParameterRange(void) const;

    /**
     * @see AbstractParametrizedCurve::operator() (double)
     */
    Eigen::Vector2d operator() (double) const;

    /**
     * @see AbstractParametrizedCurve::Derivative(double)
     */
    Eigen::Vector2d Derivative(double) const;

  private:
    /**
     * Private const field for storing the radius
     */
    const double radius;
};

#endif //PARAMETRIZEDSEMICIRCLEHPP
