///////////////////////////////////////////////////////////////////////////////
/// \file buildK.cpp
/// \brief This file provides functions to read the input-parameters, which are
///        relevant to compute the Galerkin-Matrix K corresponding to the
///        double-layer potential. The matrix is given by
///  \f[ K_{ij} = -\frac{1}{2 \pi} \int_{Ei}\int_{supp \phi_j} \frac{<y-x,n>}
///               {\vert y-x \vert^2} \phi_j(y) ds_y ds_x.  \f]
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

#include "buildK.hpp"
#include "constants.hpp"
#include "doubleLayerPotential.hpp"

void computeK(Eigen::MatrixXd& K, const Eigen::MatrixXd& coordinates,
              const Eigen::MatrixXi& elements, double eta)
{

  // resize and initialize matrix
  int nE = elements.rows();
  int nC = coordinates.rows();
  K.resize(nE, nC);
  K.setZero();
  
  // traverse the elements
  for (int j=0;j<nE;++j)
  {
    // get vertices indices and coordinates for Ei=[a,b]
    const Eigen::Vector2d& a = coordinates.row(elements(j,0));
    const Eigen::Vector2d& b = coordinates.row(elements(j,1));
   
    // traverse the elements
    for (int i=0;i<nE;++i)
    {
      // get vertices indices and coordinates for Ej=[c,d]
      int cidx = elements(i,0);
      int didx = elements(i,1);
      const Eigen::Vector2d& c = coordinates.row(cidx);
      const Eigen::Vector2d& d = coordinates.row(didx);
      
      double linetest1 = fabs( (a-c)[0]*(b-a)[1]-(a-c)[1]*(b-a)[0] );
      double linetest2 = fabs( (a-d)[0]*(b-a)[1]-(a-d)[1]*(b-a)[0] );

      if( linetest1>EPS*(a-c).norm() || linetest2>EPS*(a-d).norm() )
      {
        // compute elements' contribution
        double I0=0.0, I1=0.0;
        computeKij(&I0,&I1,eta,a,b,c,d);
        // distribute among matrix entries
        K(j,cidx) += I0-I1;
        K(j,didx) += I0+I1;
      }

    }
  }

}

