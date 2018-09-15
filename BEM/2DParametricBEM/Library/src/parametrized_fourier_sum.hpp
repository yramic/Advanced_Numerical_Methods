/**
 * \file parametrized_fourier_sum.hpp
 * \brief This file declares the class for representing Fourier Sum
 *        based parametrization. This class inherits from the
 *        Abstract base class representing parametrized curves
 * @see abstract_parametrized_curve.hpp
 *
 *  This File is a part of the 2D-Parametric BEM package
 */

#ifndef PARAMETRIZEDFOURIERSUMHPP
#define PARAMETRIZEDFOURIERSUMHPP

#include "abstract_parametrized_curve.hpp"

class ParametrizedFourierSum : public AbstractParametrizedCurve {
  public:
    using CoefficientsList = typename Eigen::MatrixXd;

    /**
     * Constructor with two coefficient list type parameters
     * that is Eigen::MatrixXd with size: 2 X N.
     * Here N is the number of sine and cosine terms in the sum
     *
     * @param cos_list Coefficient list for cosine terms
     * @param sin_list Coefficient list for sine terms
     */
    ParametrizedFourierSum(CoefficientsList cos_list,
                           CoefficientsList sin_list);

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
     * List of coefficients for the cosine terms in Fourier Sum based parametrization
     */
    const CoefficientsList cosine;

    /**
     * List of coefficients for the sine terms in Fourier Sum based parametrization
     */
    const CoefficientsList sine;
};

#endif //PARAMETRIZEDFOURIERSUMHPP
