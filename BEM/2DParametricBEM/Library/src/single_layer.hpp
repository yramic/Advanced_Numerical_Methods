/**
 * \file single_layer.hpp
 * \brief This file declares the functions to evaluate the entries of
 *        Galerkin matrices using the transformations given in section
 *        1.4.3.4 in the Lecture Notes for Advanced Numerical Methods
 *        for CSE.
 *
 * This File is a part of the 2D-Parametric BEM package
 */

#ifndef SINGLELAYERHPP
#define SINGLELAYERHPP

#include <Eigen/Dense>

#include "abstract_parametrized_curve.hpp"

namespace parametricbem2d {
  using BasisFunctionPointer = AbstractBEMSpace::BasisFunctionPointer;

  /**
   * This function is used to evaluate the Interaction Matrix for
   * the Single Layer BIO. It uses the space S_{p}^{-1} for the
   * shape functions where p is an input for the function
   *
   * @param pi Parametrization for the first panel \f$\Pi\f$.
   * @param pi_p Parametrization for the second panel \f$\Pi\f$.
   * @param p The polynomial order for the space to be used.
   * @return An Eigen::MatrixXd type Interaction Matrix (QXQ) where Q depends on p.
   */
  double ComputeIntegralAdjacent(const AbstractParametrizedCurve& pi,
                                 const AbstractParametrizedCurve& pi_p,
                                 BasisFunctionPointer bi,
                                 BasisFunctionPointer bj);

  /**
   * This function is used to evaluate the Interaction Matrix for
   * the Single Layer BIO. It uses the space S_{p}^{-1} for the
   * shape functions where p is an input for the function
   *
   * @param pi Parametrization for the first panel \f$\Pi\f$.
   * @param pi_p Parametrization for the second panel \f$\Pi\f$.
   * @param p The polynomial order for the space to be used.
   * @return An Eigen::MatrixXd type Interaction Matrix (QXQ) where Q depends on p.
   */
  double ComputeIntegralCoinciding(const AbstractParametrizedCurve& pi,
                                   const AbstractParametrizedCurve& pi_p,
                                   BasisFunctionPointer bi,
                                   BasisFunctionPointer bj);

 /**
  * This function is used to evaluate the Interaction Matrix for
  * the Single Layer BIO. It uses the space S_{p}^{-1} for the
  * shape functions where p is an input for the function
  *
  * @param pi Parametrization for the first panel \f$\Pi\f$.
  * @param pi_p Parametrization for the second panel \f$\Pi\f$.
  * @param p The polynomial order for the space to be used.
  * @return An Eigen::MatrixXd type Interaction Matrix (QXQ) where Q depends on p.
  */
  double ComputeIntegralGeneral(const AbstractParametrizedCurve& pi,
                                const AbstractParametrizedCurve& pi_p,
                                BasisFunctionPointer bi,
                                BasisFunctionPointer bj);
  /**
   * This function is used to evaluate the Interaction Matrix for
   * the Single Layer BIO. It uses the space S_{p}^{-1} for the
   * shape functions where p is an input for the function
   *
   * @param pi Parametrization for the first panel \f$\Pi\f$.
   * @param pi_p Parametrization for the second panel \f$\Pi\f$.
   * @param p The polynomial order for the space to be used.
   * @return An Eigen::MatrixXd type Interaction Matrix (QXQ) where Q depends on p.
   */
  Eigen::MatrixXd SingleLayer(const AbstractParametrizedCurve& pi,
                              const AbstractParametrizedCurve& pi_p,
                              const int& p);
} // namespace parametricbem2d

#endif //SINGLELAYERHPP
