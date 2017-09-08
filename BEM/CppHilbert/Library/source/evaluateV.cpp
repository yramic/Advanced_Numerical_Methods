///////////////////////////////////////////////////////////////////////////////
/// \file evaluateV.hpp
/// \brief This file contains the function evaluateV that evaluates the single
///        layer potential operator V tilde on any number of evaluation points
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
#include "evaluateV.hpp"


void evaluateV(Eigen::VectorXd& Vphi_x, const Eigen::MatrixXd &coordinates,
               const Eigen::MatrixXi &elements, const Eigen::VectorXd &phi,
               const Eigen::MatrixXd &x, double eta)
{
  int nX = x.rows();
  int nE = elements.rows();
  /* Initialize output vector */
  Vphi_x.resize(nX);
  Vphi_x.setZero();

  const double* gaussPoint  = getGaussPoints(GAUSS_ORDER);
  const double* gaussWeight = getGaussWeights(GAUSS_ORDER);

  /* For each boundary element Ej = [a,b], we calculate (a+b)/2
   * and (b-a)/2, because we need these vectors multiple times below.
   */
  Eigen::VectorXd lengthE(nE);
  Eigen::MatrixXd mEl(nE,2);
  Eigen::MatrixXd dEl(nE,2);
  // traverse the elements
  for (int j = 0; j < nE; ++j){
      // get vertices indices and coordinates for Ei=[a,b]
      const Eigen::Vector2d& a = coordinates.row(elements(j,0));
      const Eigen::Vector2d& b = coordinates.row(elements(j,1));
      // fill the vectors
      lengthE(j) = (b-a).norm();
      mEl.row(j)  = 0.5*(a+b);
      dEl.row(j)  = 0.5*(a-b);
  }

  for(int i = 0; i < nX; ++i){
    // save current point for readibility
    const Eigen::Vector2d& xi = x.row(i);

    for (int j = 0; j < nE; ++j){
      // get vertices indices and coordinates for Ei=[a,b]
      int aidx = elements(j,0);
      int bidx = elements(j,1);
      const Eigen::Vector2d& a = coordinates.row(aidx);
      const Eigen::Vector2d& b = coordinates.row(bidx);

      // check admissibility
      double dist_x_Ej = distancePointToSegment(xi, a, b);
      if (lengthE(j) <= eta * dist_x_Ej) {
        double sum = 0.;
        for (int k = 0; k < GAUSS_ORDER; ++k) {
            // transform point
          const Eigen::Vector2d& s = mEl.row(j) + gaussPoint[k]*dEl.row(j);
          sum += phi(j)*gaussWeight[k]*lengthE(j)*log((xi - s).squaredNorm());
        }
        Vphi_x(i) -= sum;
      }
      else {
        Vphi_x(i) -= phi(j)*lengthE(j)*slp(0, dEl.row(j), xi-mEl.row(j).transpose());
      }
    }
    Vphi_x[i] /= 8.*M_PI;
  }
}

