///////////////////////////////////////////////////////////////////////////////
/// \file buildK.hpp
/// \brief This file provides functions to read the input-parameters, which are
///        relevant to compute the Galerkin-Matrix K corresponding to the
///        double-layer potential. The matrix is given by
///  \f[ K_{ij} = -\frac{1}{2 \pi} \int_{Ei}\int_{supp \phi_j} \frac{<y-x,n>}
///               {\vert y-x \vert^2} \phi_j(y) ds_y ds_x.  \f]
///
///  This file is part of the HILBERT program package for the numerical
///  solution of the Laplace equation with mixed boundary conditions by use of
///  BEM in 2D.
///
///  C++ adaptation for ANCSE17 of HILBERT V3.1 TUWien 2009-2013
///////////////////////////////////////////////////////////////////////////////
#ifndef _BUILDK_HPP_
#define _BUILDK_HPP_

#include <Eigen/Dense>

/**
 *  This function generates the panels Ei and Ej. In case of Ej and Ei being
 *  segments of the same line the entry K_ij vanishes. Otherwise the function
 *  computeKij() is called for the computation of the corresponding matrix
 *  entry.
 *
 *  @param[out] K  nE x nC matrix. nC and nE are the number of coordinates and
 *                 elements, respectively.
 *  @param[in] coordinates  nC x 2 matrix containing the coordinates of the
 *                          vertices of the boundary mesh.
 *  @param[in] elements  nE x 2 matrix containing the indices of the vertices
 *                       corresponding to each element of the boundary mesh.
 *  @param[in] eta  Admissibility constant. It is greater or equal than 0.
 */
void computeK(Eigen::MatrixXd& K, const Eigen::MatrixXd& coordinates,
              const Eigen::MatrixXi& elements, double eta);

#endif

