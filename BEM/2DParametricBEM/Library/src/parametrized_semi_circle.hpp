/*
 * This hpp file declares the class representing parametrization
 * of a semi circle centered at (0,0), in the positive x-half
 * with the given radius. This class is inherited from the
 * Abstract base class representing parametrized curves
 * (See abstract_parametrized_curve.hpp)
 */

#ifndef PARAMETRIZEDSEMICIRCLEHPP
#define PARAMETRIZEDSEMICIRCLEHPP

#include <utility>

#include <Eigen/Dense>

#include "abstract_parametrized_curve.hpp"

class ParametrizedSemiCircle : public AbstractParametrizedCurve {
  public:
    /* Constructor with specified radius; default value = 1.*/
    ParametrizedSemiCircle(double = 1.);

    std::pair<double,double> ParameterRange(void) const;

    Eigen::Vector2d operator() (double) const;

    Eigen::Vector2d Derivative(double) const;

  private:
    /* Private const field for radius*/
    const double radius;
};

#endif //PARAMETRIZEDSEMICIRCLEHPP
