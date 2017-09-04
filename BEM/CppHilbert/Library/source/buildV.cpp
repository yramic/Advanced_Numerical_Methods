///////////////////////////////////////////////////////////////////////////////
/// \file buildV.cpp
/// \brief This file provides functions for the computation of the single-layer
///  potential matrix V defined by
///  \f[ V_{jk} = < V \chi_j , \chi_k >  \f]
///  with characteristic functions \f$\chi_j\f$ corresponding to a line segment
///  Tj.
///  The single-layer potential is defined by
///  \f[ V \chi_j = -\frac{1}{2\pi} \int_{Ej} \log{ \vert x-y \vert } ds_y. \f]
///
///  This file contains only the implementation. For extensive documentation
///  consult the corresponding header-file.
///
///  This file is part of the HILBERT program package for the numerical
///  solution of the Laplace equation with mixed boundary conditions by use of
///  BEM in 2D.
///
///  C++ adaptation for ANCSE17 of HILBERT V3.1 TUWien 2009-2013
///////////////////////////////////////////////////////////////////////////////
#include <cmath>

#include "buildV.hpp"
#include "constants.hpp"
#include "singleLayerPotential.hpp"


//------------------------------------------------------------------------------
void computeV(Eigen::MatrixXd& V, const Eigen::MatrixXd& coordinates,
	      const Eigen::MatrixXi& elements, double eta)
{

  assert(eta >= 0);

  // resize matrix
  int nE = elements.rows();
  V.resize(nE,nE);
  // traverse the elements
  for (int i=0; i<nE; ++i)
  {
    // get vertices indices and coordinates for Ei=[a,b]
    const Eigen::Vector2d& a = coordinates.row(elements(i,0));
    const Eigen::Vector2d& b = coordinates.row(elements(i,1));

    // traverse the elements
    for (int j=i; j<nE; ++j)
    {
      // get vertices indices and coordinates for Ej=[c,d]
      const Eigen::Vector2d& c = coordinates.row(elements(j,0));
      const Eigen::Vector2d& d = coordinates.row(elements(j,1));

      // compute elements' contribution
      double tmp = computeVij(a, b, c, d, eta);
      // distribute it among the matrix entries
      V(i,j) = tmp;
      V(j,i) = V(i,j);
    }
  }
  
}

