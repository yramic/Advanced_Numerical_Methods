///////////////////////////////////////////////////////////////////////////////
/// \file evaluateV.hpp
/// \brief This file contains the function evaluateV that evaluates the single
///        layer potential operator V tilde on any number of evaluation points
///        within the domain \f$Omega\f$.
///
///  This file is part of the HILBERT program package for the numerical
///  solution of the Laplace equation with mixed boundary conditions by use of
///  BEM in 2D.
///
///  C++ adaptation for ANCSE17 of HILBERT V3.1 TUWien 2009-2013
///////////////////////////////////////////////////////////////////////////////
#ifndef _EVALUATEV_HPP_GUARD_
#define _EVALUATEV_HPP_GUARD_

#include "singleLayerPotential.hpp"
#include "geometry.hpp"

extern "C" {
#include "gaussQuadrature.h"
}

/**
 *  evaluateV evaluates the single layer potential operator V tilde on a number
 *  of evaluation points.
 *
 *  @param[out] Vphi_x vector
 *  @param[in] coordinates  (nC x 2) matrix containing the coordinates of the
 *                          vertices of the boundary mesh.
 *  @param[in] elements  (nE x 2) matrix containing the indices of the vertices
 *                       corresponding to each element of the boundary mesh.
 *  @param[in] phi nC vector such that phi_i=phi(z_i).
 *  @param[in] x  (nX x 2) matrix that contains the evaluation points.
 *  @param[in] eta  Admissibility constant. It is greater or equal than 0.
 */
void evaluateV(Eigen::VectorXd& Vphi_x, const Eigen::MatrixXd& coordinates,
               const Eigen::MatrixXi& elements, const Eigen::VectorXd& phi,
               const Eigen::MatrixXd& x, double eta);

#endif

