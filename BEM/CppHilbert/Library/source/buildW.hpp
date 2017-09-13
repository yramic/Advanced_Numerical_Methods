///////////////////////////////////////////////////////////////////////////////
/// \file buildW.hpp
/// \brief This file provides functions for the calculation of the
///  hypersingular integral operator matrix W that is defined by
///  \f[ W_{jk} = < W \phi_j, \phi_k >  \f]
///  with the hat functions \f$\phi_i\f$ corresponding to a piece-wise affine
///  function that is linear on any element of the boundary mesh and that is 1
///  on the i-th node z_i and 0 on any other node.
///
///  The matrix W is in fact calculated using the following relation between W
///  and the simple-layer potential V:
///  \f[ <Wu, v> = <Vu', v'>   \f]
///  In this formula, u' and v' denote the arc-length derivative of u and v,
///  respectively.
///
///  This file is part of the HILBERT program package for the numerical
///  solution of the Laplace equation with mixed boundary conditions by use of
///  BEM in 2D.
///
///  C++ adaptation for ANCSE17 of HILBERT V3.1 TUWien 2009-2013
///////////////////////////////////////////////////////////////////////////////
#ifndef _BUILDW_HPP_GUARD_
#define _BUILDW_HPP_GUARD_

#include <Eigen/Dense>


/**
 *  This function assembles the matrix W by calling computeWij() for each pair
 *  of boundary elements (Ei, Ej) and storing the result in W. It uses the
 *  symmetry of W to reduce the build time.
 *
 *  @param[out] W  (nC x nC) matrix. nC and nE are the number of coordinates and
 *                 elements, respectively.
 *  @param[in] coordinates  (nC x 2) matrix containing the coordinates of the
 *                          vertices of the boundary mesh.
 *  @param[in] elements  (nE x 2) matrix containing the indices of the vertices
 *                       corresponding to each element of the boundary mesh.
 *  @param[in] eta  Admissibility constant. It is greater or equal than 0.
 */
void computeW(Eigen::MatrixXd& W, const Eigen::MatrixXd& coordinates,
	      const Eigen::MatrixXi& elements, double eta);

#endif

