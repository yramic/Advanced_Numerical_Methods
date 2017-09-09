///////////////////////////////////////////////////////////////////////////////
/// \file buildV.hpp
/// \brief This file provides functions for the computation of the single-layer
///  potential matrix V defined by
///  \f[ V_{jk} = < V \chi_j , \chi_k >  \f]
///  with characteristic functions \f$\chi_j\f$ corresponding to a line segment
///  Tj.
///  The single-layer potential is defined by
///  \f[ V \chi_j = -\frac{1}{2\pi} \int_{Ej} \log{ \vert x-y \vert } ds_y. \f]
///
///  This file is part of the HILBERT program package for the numerical
///  solution of the Laplace equation with mixed boundary conditions by use of
///  BEM in 2D.
///
///  C++ adaptation for ANCSE17 of HILBERT V3.1 TUWien 2009-2013
///////////////////////////////////////////////////////////////////////////////
#ifndef _BUILDV_HPP_GUARD_
#define _BUILDV_HPP_GUARD_

#include <Eigen/Dense>


/**
 *  This function assembles the matrix V by calling computeVij() for each pair
 *  of boundary elements (Ei, Ej) and storing the result in V. It uses the
 *  symmetry of V to reduce the build time.
 *
 *  @param[out] V  (nE x nE) matrix. nC and nE are the number of coordinates and
 *                 elements, respectively.
 *  @param[in] coordinates  (nC x 2) matrix containing the coordinates of the
 *                          vertices of the boundary mesh.
 *  @param[in] elements  (nE x 2) matrix containing the indices of the vertices
 *                       corresponding to each element of the boundary mesh.
 *  @param[in] eta  Admissibility constant. It is greater or equal than 0.
 */
void computeV(Eigen::MatrixXd& V, const Eigen::MatrixXd& coordinates,
	      const Eigen::MatrixXi& elements, double eta);

#endif

