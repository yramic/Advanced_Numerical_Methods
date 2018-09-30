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
#include "abstract_bem_space.hpp"
#include "parametrized_mesh.hpp"

namespace parametricbem2d {
  using BasisFunctionPointer = AbstractBEMSpace::BasisFunctionPointer;

  /**
   * This function is used to evaluate a single entry of the
   * Interaction Matrix for the pair of panels \f$\Pi\f$ and \f$\Pi\f$'
   * for the bilinear form induced by the Single Layer BIO. It implements
   * the case where the panels are adjacent. This function calculates the
   * matrix entry by using a local arclength parametrization and the
   * transformations mentioned in eq. 1.4.165 to eq. 1.4.172
   *
   * @param pi Parametrization for the first panel \f$\Pi\f$.
   * @param pi_p Parametrization for the second panel \f$\Pi\f$'.
   * @param bi The basis function associated with panel \f$\Pi\f$.
   * @param bj The basis function associated with panel \f$\Pi\f$'.
   * @return The (i,j)^{th} entry of the interaction matrix
   */
  double ComputeIntegralAdjacent(const AbstractParametrizedCurve& pi,
                                 const AbstractParametrizedCurve& pi_p,
                                 BasisFunctionPointer bi,
                                 BasisFunctionPointer bj);

  /**
   * This function is used to evaluate a single entry of the
   * Interaction Matrix for the pair of panels \f$\Pi\f$ and \f$\Pi\f$'
   * for the bilinear form induced by the Single Layer BIO. It implements
   * the case where the panels are coinciding. The value is calculated
   * using the transformations mentioned in eq. 1.4.153 to eq. 1.4.163
   *
   * @param pi Parametrization for the first panel \f$\Pi\f$.
   * @param pi_p Parametrization for the second panel \f$\Pi\f$'.
   * @param bi The basis function associated with panel \f$\Pi\f$.
   * @param bj The basis function associated with panel \f$\Pi\f$'.
   * @return The (i,j)^{th} entry of the interaction matrix
   */
  double ComputeIntegralCoinciding(const AbstractParametrizedCurve& pi,
                                   const AbstractParametrizedCurve& pi_p,
                                   BasisFunctionPointer bi,
                                   BasisFunctionPointer bj);

 /**
  * This function is used to evaluate a single entry of the
  * Interaction Matrix for the pair of panels \f$\Pi\f$ and \f$\Pi\f$'
  * for the bilinear form induced by the Single Layer BIO. It implements
  * the general case where panels are neither coinciding nor adjacent.
  * The value is calculated using standard Gaussian Quadrature.
  *
  * @param pi Parametrization for the first panel \f$\Pi\f$.
  * @param pi_p Parametrization for the second panel \f$\Pi\f$'.
  * @param bi The basis function associated with panel \f$\Pi\f$.
  * @param bj The basis function associated with panel \f$\Pi\f$'.
  * @return The (i,j)^{th} entry of the interaction matrix
  */
  double ComputeIntegralGeneral(const AbstractParametrizedCurve& pi,
                                const AbstractParametrizedCurve& pi_p,
                                BasisFunctionPointer bi,
                                BasisFunctionPointer bj);
  /**
   * This function is used to evaluate the Interaction Matrix for
   * the pair of panels \f$\Pi\f$ and \f$\Pi\f$' for the bilinear
   * form induced by the Single Layer BIO given by the formula :
   * \f$I_{ij}\f$ = \f$-\frac{1}{2\pi} \int_{\Pi} \int_{\Pi '} \log{ \vert x-y \vert } b^{j}_{\Pi '}(y) b^{i}_{\Pi}(x) ds(y) ds(x) \f$
   * where \f$b^{j}_{\Pi '}\f$ & \f$b^{i}_{\Pi}\f$ are local shape functions
   * associated with panels \f$\Pi\f$ ' and \f$\Pi\f$ respectively. \f$I\f$, the
   * interaction matrix is of size QXQ where Q is the number of local shape
   * functions for the used BEM space. The computation of the entries
   * are based on cases and delegated to these functions accordingly:
   * ComputeIntegralGeneral
   * ComputeIntegralAdjacent
   * ComputeIntegralCoinciding
   *
   * @param pi Parametrization for the first panel \f$\Pi\f$.
   * @param pi_p Parametrization for the second panel \f$\Pi\f$'.
   * @param space The BEM space to be used for calculations
   * @return An Eigen::MatrixXd type Interaction Matrix (QXQ)
   *         where Q depends on the BEM space employed.
   */
  Eigen::MatrixXd SingleLayer(const AbstractParametrizedCurve& pi,
                              const AbstractParametrizedCurve& pi_p,
                              const AbstractBEMSpace& space);

  /**
   * This function is used to evaluate the Interaction Matrix for
   * the pair of panels \f$\Pi\f$ and \f$\Pi\f$' for the bilinear
   * form induced by the Single Layer BIO given by the formula :
   * \f$I_{ij}\f$ = \f$-\frac{1}{2\pi} \int_{\Pi} \int_{\Pi '} \log{ \vert x-y \vert } b^{j}_{\Pi '}(y) b^{i}_{\Pi}(x) ds(y) ds(x) \f$
   * where \f$b^{j}_{\Pi '}\f$ & \f$b^{i}_{\Pi}\f$ are local shape functions
   * associated with panels \f$\Pi\f$ ' and \f$\Pi\f$ respectively. \f$I\f$, the
   * interaction matrix is of size QXQ where Q is the number of local shape
   * functions for the used BEM space. The computation of the entries
   * are based on cases and delegated to these functions accordingly:
   * ComputeIntegralGeneral
   * ComputeIntegralAdjacent
   * ComputeIntegralCoinciding
   *
   * @param pi Parametrization for the first panel \f$\Pi\f$.
   * @param pi_p Parametrization for the second panel \f$\Pi\f$'.
   * @param space The BEM space to be used for calculations
   * @return An Eigen::MatrixXd type Interaction Matrix (QXQ)
   *         where Q depends on the BEM space employed.
   */
  Eigen::MatrixXd SingleLayerMatrix(const ParametrizedMesh mesh,
                                    const AbstractBEMSpace& space);
} // namespace parametricbem2d

#endif //SINGLELAYERHPP
