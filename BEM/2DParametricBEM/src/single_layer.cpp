/**
 * \file single_layer.cpp
 * \brief This file defines the functions to evaluate the entries of
 *        Galerkin matrices using the transformations given in section
 *        1.4.3.4 in the Lecture Notes for Advanced Numerical Methods
 *        for CSE.
 *
 * This File is a part of the 2D-Parametric BEM package
 */

#include "single_layer.hpp"

#include <math.h>
#include <vector>

#include <Eigen/Dense>
#include "abstract_parametrized_curve.hpp"
#include "abstract_bem_space.hpp"
#include "discontinuous_space.hpp"
extern "C" {
#include "gaussQuadrature.h"
}
#include "logweight_quadrature.hpp"
//#include "integration.h"

namespace parametricbem2d {

  Eigen::MatrixXd SingleLayer(const AbstractParametrizedCurve& pi,
                              const AbstractParametrizedCurve& pi_p,
                              const AbstractBEMSpace& space) {
    int Q = space.getQ(); // The number of Reference Shape Functions in space
    Eigen::MatrixXd interaction_matrix(Q,Q);
    // Vector containing the Reference Shape Functions
    std::vector<BasisFunctionPointer> bases = space.getShapeFunctions();
    unsigned int N = 32; // No. of points for Gauss Quadrature

    for (int i=0 ; i<Q ; ++i) {
      for (int j=0; j<Q ; ++j) {
        if (&pi==&pi_p) // Same Panels
          interaction_matrix(i,j) = ComputeIntegralCoinciding(pi,
                                                              pi_p,
                                                              bases[i],
                                                              bases[j]);
        else if ( pi(1)==pi_p(-1) || pi(-1)==pi_p(1)) // Adjacent Panels
          interaction_matrix(i,j) = ComputeIntegralAdjacent(pi,
                                                            pi_p,
                                                            bases[i],
                                                            bases[j]);
        else
          interaction_matrix(i,j) = ComputeIntegralGeneral(pi,
                                                           pi_p,
                                                           bases[i],
                                                           bases[j]);
      }
    }
    return interaction_matrix;
  }

  double ComputeIntegralCoinciding(const AbstractParametrizedCurve& pi,
                                   const AbstractParametrizedCurve& pi_p,
                                   BasisFunctionPointer bi,
                                   BasisFunctionPointer bj) {
    unsigned int N = 32;
    // Lambda expression for functions F and G
    auto F = [&] (double t) {
      return bj(t)*pi_p.Derivative(t).norm();
    };

    auto G = [&] (double s) {
      return bi(s)*pi.Derivative(s).norm();
    };

    // Lambda expression for the 1st integrand in eq. 1.4.159
    auto integrand1 = [&] (double s, double t) {
      double tol = 1e-5;
      double s_st;
      if (fabs(s-t)>tol)
        s_st = (pi(s)-pi_p(t)).squaredNorm()/(s-t)/(s-t);
      else
        s_st = pi.Derivative(t).squaredNorm();
      return 0.5*log(s_st)*F(t)*G(s);
    };

    double i1=0.,i2=0.; // The two integrals in eq. 1.4.159

    // Getting quadrature weights and points
    const double* weights = getGaussWeights(N);
    const double* points = getGaussPoints(N);

    // Double sum for double 1st integral in eq. 1.4.159
    for ( unsigned int i = 0; i<N ; ++i ) {
      for ( unsigned int j = 0; j<N ; ++j ) {
        i1 += weights[i]*weights[j]*integrand1(points[i],points[j]);
      }
    }

    // Lambda expression for 2nd integrand (1.4.159) in transformed
    // coordinates (1.4.163)
    auto integrand2 = [&] (double w, double z) {
      return F(0.5*(w-z))*G(0.5*(w+z))+F(0.5*(w+z))*G(0.5*(w-z));
    };

    // Getting log weighted quadrature nodes and weights
    QuadRule logweightQR = getLogWeightQR(2.,N);

    // Double sum for double 2nd integral (1.4.163)
    for ( unsigned int i = 0; i<N ; ++i ) {
      double z = logweightQR.x(i);
      double inner = 0.;
      for ( unsigned int j = 0; j<N ; ++j ) {
        double w = points[j] * (2-z);
        inner += weights[j] * integrand2(w,z);
      }
      inner *= (2-z);
      i2 += logweightQR.w(i) * inner;
    }
    return -1./(2.*M_PI)*(i1+0.5*i2);
  }

  double ComputeIntegralAdjacent(const AbstractParametrizedCurve& pi,
                                 const AbstractParametrizedCurve& pi_p,
                                 BasisFunctionPointer bi,
                                 BasisFunctionPointer bj) {
    unsigned int N = 32;
    // Getting the panel lengths for local arclength parametrization
    double length_pi = 1.; // Panel length for local arclength parametrization
    double length_pi_p = 1.; // Panel length for local arclength parametrization
    // swap ensures when transforming the parametrizations from [-1,1]->\Pi to
    // local arclength parametrizations [0,|\Pi|] -> \Pi, the common point
    // between the panels corresponds to the parameter 0 in arclength
    // parametrizations
    bool swap = (pi(1)!=pi_p(-1));
    auto F = [&] (double t_pr) { // Function associated with panel pi_p
      double t = swap ? 1-2*t_pr/length_pi_p : 2*t_pr/length_pi_p-1;
      return bj(t)*pi_p.Derivative(t).norm();
    };
    auto G = [&] (double s_pr) { // Function associated with panel pi
      double s = swap ? 2*s_pr/length_pi-1 : 1-2*s_pr/length_pi;
      return bi(s)*pi.Derivative(s).norm();
    };
    auto D_r_phi = [&] (double r,double phi) { // eq. 1.4.172
      double s_pr = r * cos(phi); double s = swap ? 2*s_pr/length_pi-1 : 1-2*s_pr/length_pi;
      double t_pr = r * sin(phi); double t = swap ? 1-2*t_pr/length_pi_p : 2*t_pr/length_pi_p-1;
      if (r!=0)
        return (pi(s)-pi_p(t)).squaredNorm()/r/r;
      else
        return 1-sin(2*phi)*pi.Derivative(s).dot(pi_p.Derivative(t));
    };

    double alpha = atan(length_pi_p/length_pi);
    // Getting Gauss Quadrature weights and nodes
    const double* weights = getGaussWeights(N);
    const double* points = getGaussPoints(N);

    // The two integrals in eq. 1.4.172 have to be further split into two parts
    // part 1 is where phi goes from 0 to alpha
    // part 2 is where phi goes from alpha to pi/2

    // i_IJ -> Integral I, part J
    double i11 = 0.,i21 = 0.,i12 = 0., i22 = 0.;
    // part 1
    for (unsigned int i = 0; i<N ; ++i) {
      // Transforming gauss quadrature abcissa into phi (0 to alpha)
      double phi = alpha/2*(1+points[i]);
      // Computing log weighted inner integral with fixed phi
      double inner1=0.; // inner integral for double integral 1
      double inner2=0.; // inner integral for double integral 2
      double rmax = length_pi/cos(phi); // Upper limit for inner 'r' integral
      for (unsigned int j = 0; j<N ; ++j) {
        // Getting Quadrature weights and nodes for log weight
        QuadRule logweightQR = getLogWeightQR(rmax,N);
        double r = logweightQR.x(j); //for inner2 log weighted quadrature
        inner2 += logweightQR.w(j) * r * F(r*sin(phi)) * G(r*cos(phi));
        r = rmax/2*(1+points[j]); //for inner1 standard gauss quadrature
        inner1 += weights[j] * r * log(D_r_phi(r,phi)) * F(r*sin(phi)) * G(r*cos(phi));
      }
      inner1 *= rmax/2;
      i11 += weights[i] * inner1 * alpha/2;
      i21 += weights[i] * inner2 * alpha/2;
    }

    // part 2 (phi from alpha to pi/2)
    for (unsigned int i = 0; i<N ; ++i) {
      // Transforming gauss quadrature abcissa into phi (alpha to pi/2)
      double phi = points[i]*(M_PI/2.-alpha)/2.+(M_PI/2.+alpha)/2.;
      // Computing log weighted inner integral with fixed phi
      double inner1=0.; // inner integral for double integral 1
      double inner2=0.; // inner integral for double integral 2
      double rmax = length_pi_p/sin(phi); // Upper limit for inner 'r' integral
      for (unsigned int j = 0; j<N ; ++j) {
        // Getting Quadrature weights and nodes for log weight
        QuadRule logweightQR = getLogWeightQR(rmax,N);
        double r = logweightQR.x(j); //for inner2 log weighted quadrature
        inner2 += logweightQR.w(j) * r * F(r*sin(phi)) * G(r*cos(phi));
        r = rmax/2*(1+points[j]); //for inner1 standard gauss quadrature
        inner1 += weights[j] * r * log(D_r_phi(r,phi)) * F(r*sin(phi)) * G(r*cos(phi));
      }
      inner1 *= rmax/2;
      i12 += weights[i] * inner1 *(M_PI/2.-alpha)/2.;
      i22 += weights[i] * inner2 *(M_PI/2.-alpha)/2.;
    }
    double integral = 0.5*(i11+i12) + (i21+i22);
    integral *= 4/length_pi/length_pi_p;
    return -1/(2*M_PI)*integral;
  }

  double ComputeIntegralGeneral(const AbstractParametrizedCurve& pi,
                                const AbstractParametrizedCurve& pi_p,
                                BasisFunctionPointer bi,
                                BasisFunctionPointer bj) {
    // Lambda expression for functions F and G
    auto F = [&] (double t) {
      return bj(t)*pi_p.Derivative(t).norm();
    };

    auto G = [&] (double s) {
      return bi(s)*pi.Derivative(s).norm();
    };

    /*alglib::ae_int_t N = 100;
    alglib::ae_int_t info;
    alglib::real_1d_array points,weights;
    alglib::gqgenerategausslegendre(N,info,points,weights);*/

    double integral = 0.;
    unsigned int N = 32;
    // Getting quadrature weights and points
    const double* weights = getGaussWeights(N);
    const double* points = getGaussPoints(N);

    // Double sum for double integral
    for ( unsigned int i = 0; i<N ; ++i ) {
      for ( unsigned int j = 0; j<N ; ++j ) {
        double s = points[i];
        double t = points[j];
        integral += weights[i]*weights[j]*log((pi(s)-pi_p(t)).norm())*F(t)*G(s);
      }
    }
    return -1/(2*M_PI)*integral;
  }


} // namespace parametricbem2d
