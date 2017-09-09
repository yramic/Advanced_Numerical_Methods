///////////////////////////////////////////////////////////////////////////////
/// \file evaluateK.hpp
/// \brief This file contains the function evaluateK that evaluates the double
///        layer potential operator K tilde on any number of evaluation points
///        within the domain \f$Omega\f$.
///
///  This file is part of the HILBERT program package for the numerical
///  solution of the Laplace equation with mixed boundary conditions by use of
///  BEM in 2D.
///
///  C++ adaptation for ANCSE17 of HILBERT V3.1 TUWien 2009-2013
///////////////////////////////////////////////////////////////////////////////
#ifndef _EVALUATEK_HPP_GUARD_
#define _EVALUATEK_HPP_GUARD_

#include "geometry.hpp"
#include "doubleLayerPotential.hpp"
extern "C" {
#include "gaussQuadrature.h"
}


/**
 *  evaluateK evaluates the double layer potential operator K tilde on a number
 *  of evaluation points.
 *
 *  @param[out] Kgx vector
 *  @param[in] coordinates  (nC x 2) matrix containing the coordinates of the
 *                          vertices of the boundary mesh.
 *  @param[in] elements  (nE x 2) matrix containing the indices of the vertices
 *                       corresponding to each element of the boundary mesh.
 *  @param[in] gh nC vector such that gh_i=g(z_i).
 *  @param[in] x  (nX x 2) matrix that contains the evaluation points.
 *  @param[in] eta  Admissibility constant. It is greater or equal than 0.
 */
void evaluateK(Eigen::VectorXd& Kgx, const Eigen::MatrixXd& coordinates,
               const Eigen::MatrixXi& elements, const Eigen::VectorXd& gh,
               const Eigen::MatrixXd& x, double eta);

#endif
