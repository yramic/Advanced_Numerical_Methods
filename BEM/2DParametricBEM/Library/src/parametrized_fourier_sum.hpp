/*
 * This hpp file declares the class representing Fourier
 * Sum based parametrization which is inherited from the
 * Abstract base class representing parametrized curves
 * (See abstract_parametrized_curve.hpp)
 */

#ifndef PARAMETRIZEDFOURIERSUMHPP
#define PARAMETRIZEDFOURIERSUMHPP

#include "abstract_parametrized_curve.hpp"

class ParametrizedFourierSum : public AbstractParametrizedCurve {
  public:
    using CoefficientsList = typename Eigen::MatrixXd;

    /* Constructor with two coefficient list parameters
    * in the form of Eigen::MatrixXd (Size: 2 X N)
    * Here N is the number of sine and cosine terms in the sum
    */
    ParametrizedFourierSum(CoefficientsList,
                           CoefficientsList);

    std::pair<double,double> ParameterRange(void) const;

    Eigen::Vector2d operator() (double) const;

    Eigen::Vector2d Derivative(double) const;

  private:
    /* List of coefficients for the cosine terms in Fourier Sum based parametrization*/
    const CoefficientsList cosine;

    /* List of coefficients for the sine terms in Fourier Sum based parametrization*/
    const CoefficientsList sine;
};

#endif //PARAMETRIZEDFOURIERSUMHPP
