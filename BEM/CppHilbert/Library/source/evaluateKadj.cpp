///////////////////////////////////////////////////////////////////////////////
/// \file evaluateKadj.cpp
/// \brief This file contains the function evaluateKadj that evaluates the
///        adjoint double layer potential operator K* tilde on any number of
///        evaluation points on the boundary \f$Gamma\f$.
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
#include "doubleLayerPotential.hpp"
#include "evaluateKadj.hpp"
#include "geometry.hpp"

extern "C" {
#include "gaussQuadrature.h"
}



void evaluateKadj(Eigen::VectorXd &Kx, const Eigen::MatrixXd &coordinates,
                  const Eigen::MatrixXi &elements, const Eigen::VectorXd &gh,
                  const Eigen::MatrixXd &x, const Eigen::MatrixXd &n_x, double eta)
{

  int n_gauss_points=GAUSS_ORDER;
  const double* q_points  = getGaussPoints(GAUSS_ORDER);
  const double* q_weights = getGaussWeights(GAUSS_ORDER);

  int nX = x.rows();
  int nE = elements.rows();
  /* Initialize output vector */
  Kx.resize(nX);
  Kx.setZero();

  // traverse evaluation points
  for (int j=0;j<nX; ++j){
    // save current point for readibility
    const Eigen::Vector2d& xj = x.row(j);

    //traverse elements
    for (int i=0;i<nE;++i){
      // get vertices indices and coordinates for Ei=[a,b]
      int aidx = elements(j,0);
      int bidx = elements(j,1);
      const Eigen::Vector2d& a = coordinates.row(aidx);
      const Eigen::Vector2d& b = coordinates.row(bidx);

      // middle point and auxiliary vectors
      const Eigen::Vector2d& mp = 0.5*(a+b);
      const Eigen::Vector2d& vc = (b-a);
      const Eigen::Vector2d& u  = 0.5*vc;
      const Eigen::Vector2d& v  = xj - mp;
      // auxiliary constants
      double u_sqnorm = u.squaredNorm();
      double v_sqnorm = v.squaredNorm();
      double sqrdelta = sqrt(u_sqnorm*v_sqnorm - std::pow(u.dot(v),2));

      /*Determine the normal vector*/
      Eigen::Vector2d n_vc;
      n_vc<< (0.5 * vc[1])/sqrt(u_sqnorm), -(0.5 * vc[0])/sqrt(u_sqnorm);

      double K_0 = 0;

      /*Check if element is admissible*/
      double dist_x_Ej = distancePointToSegment(xj, a, b);
      if (sqrt(2*u_sqnorm) > eta*dist_x_Ej){ /* x and Ej are inadmissable. */
        /*Compute the integrals analytically*/

        /* It is not safe to call dlp if u, v have the same length and
         * are linearly dependent. But fortunately, this is only the case
         * if the whole integral is 0.
         */

        if (sqrdelta > EPS * sqrt(u_sqnorm * v_sqnorm)){
          double g_0m1 = dlp(0,u,v); /*Compute g_0^(-1)*/
          double g_1m1 = dlp(1,u,v); /*Compute g_1^(-1)*/

          /*Compute K_0*/
          K_0 = (n_x.row(j).dot(u)* g_1m1 + n_x.row(j).dot(v) * g_0m1);
          K_0 /= M_PI;
        }
        else{
          K_0 = 0.;
        }
      }

      else{
        /*Compute the integrals via quadrature rule*/
        K_0=0;
        for (int k=0; k<n_gauss_points; ++k){
          double contrib = n_x.row(j).dot(u*q_points[k]+v)/
                          ((u_sqnorm*q_points[k]*q_points[k]
                            + 2.*u.dot(v)*q_points[k] + v_sqnorm)*M_PI);
          K_0 += q_weights[k] * contrib;
        }
      }

      Kx[j] -= 0.5 * gh[i] * sqrt(u_sqnorm) * K_0;
    }
  }
}

