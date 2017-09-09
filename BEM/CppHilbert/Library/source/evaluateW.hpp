///////////////////////////////////////////////////////////////////////////////
/// \file evaluateW.hpp
/// \brief This file contains the function evaluateW that evaluates the
///        hypersingular integral operator W on any number of evaluation points
///        on the boundary \f$Gamma\f$.
///
///  This file is part of the HILBERT program package for the numerical
///  solution of the Laplace equation with mixed boundary conditions by use of
///  BEM in 2D.
///
///  C++ adaptation for ANCSE17 of HILBERT V3.1 TUWien 2009-2013
///////////////////////////////////////////////////////////////////////////////
#ifndef _EVALUATEW_HPP_GUARD_
#define _EVALUATEW_HPP_GUARD_


/**
 *  evaluateW evaluates the hypersingular integral operator W on a number of
 *  evaluation points.
 *
 *  @param[out] Wx vector
 *  @param[in] coordinates  (nC x 2) matrix containing the coordinates of the
 *                          vertices of the boundary mesh.
 *  @param[in] elements  (nE x 2) matrix containing the indices of the vertices
 *                       corresponding to each element of the boundary mesh.
 *  @param[in] gh nC vector such that gh_i=g(z_i).
 *  @param[in] x  (nX x 2) matrix that contains the evaluation points.
 *  @param[in] n_x (nX x 2) matrix that contains the normal vectors.
 *  @param[in] eta  Admissibility constant. It is greater or equal than 0.
 */
void evaluateW(Eigen::VectorXd& Wx, const Eigen::MatrixXd& coordinates,
               const Eigen::MatrixXi& elements, const Eigen::VectorXd& gh,
               const Eigen::MatrixXd& x,const Eigen::MatrixXd& n_x, double eta);

#endif


