///////////////////////////////////////////////////////////////////////////////
/// \file evaluateK.cpp
/// \brief This file contains the function evaluateK that evaluates the double
///        layer potential operator K tilde on any number of evaluation points
///        within the domain \f$Omega\f$.
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

#include "constants.hpp"
#include "evaluateK.hpp"


void evaluateK(Eigen::VectorXd& Kgx, const Eigen::MatrixXd& coordinates,
               const Eigen::MatrixXi& elements, const Eigen::VectorXd &gh,
               const Eigen::MatrixXd &x, double eta)
{
  int nX = x.rows();
  int nE = elements.rows();
  // Initialize output vector
  Kgx.resize(nX);
  Kgx.setZero();

  // Get quadrature points and weights
  const double* qp  = getGaussPoints(GAUSS_ORDER);
  const double* qw = getGaussWeights(GAUSS_ORDER);

  // Auxiliary vector containing elements info.
  Eigen::VectorXd lengthE(nE);
  Eigen::MatrixXd mEl(nE,2);
  Eigen::MatrixXd dEl(nE,2);
  Eigen::MatrixXd nEl(nE,2);

  // Traverse the elements
  for (int j = 0; j < nE; ++j){
    // Get vertices indices and coordinates for Ei=[a,b]
    const Eigen::Vector2d& a = coordinates.row(elements(j,0));
    const Eigen::Vector2d& b = coordinates.row(elements(j,1));
    // Fill the vectors
    lengthE(j) = (b-a).norm();
    mEl.row(j)  = 0.5*(a+b);
    dEl.row(j)  = 0.5*(b-a);
    nEl.row(j)  = unitNormal(a,b) ;
  }

  // Traverse the evaluation points
  for (int i = 0; i < nX; ++i){
      // Save current point for readibility
      const Eigen::Vector2d& xi = x.row(i);

    // Traverse elements
    for (int j = 0; j < nE; ++j){
      // Get vertices indices and coordinates for Ei=[a,b]
      int aidx = elements(j,0);
      int bidx = elements(j,1);
      const Eigen::Vector2d& a = coordinates.row(aidx);
      const Eigen::Vector2d& b = coordinates.row(bidx);

      // Check admissibility
      double dist_xi_Ej = distancePointToSegment(xi, a, b);
      if (lengthE(j) <= dist_xi_Ej * eta){
        // Integrate
        double sum = 0.;
        for (int k = 0; k < GAUSS_ORDER; ++k){
          // Transform point
          const Eigen::Vector2d& s = mEl.row(j) + qp[k]*dEl.row(j);
          // Add contribution
          sum += qw[k]*(s-xi).dot(nEl.row(j))
                  *( gh[aidx] + ((1+qp[k])*(gh[bidx]-gh[aidx])/2) )
                  / (s-xi).squaredNorm();
        }

        Kgx(i) -= lengthE(j) * sum;
      }
      else{
        // Compute integral analitically
        const Eigen::Vector2d& aux = mEl.row(j)-xi.transpose();
        if( fabs(aux.dot(nEl.row(j))) > EPS ){
          double commonFactor = aux.dot(nEl.row(j)) * lengthE(j);
          double prod1 = commonFactor * dlp(0, dEl.row(j), aux ) / 2.;
          double prod2 = commonFactor * dlp(1, dEl.row(j), aux ) / 2.;
          Kgx(i) -= ( (gh(aidx) + gh(bidx))* prod1 + (gh(bidx) - gh(aidx))* prod2);
        }
      }

    } // end elements' for loop

    Kgx(i) /= 4. * M_PI;
  }

}

