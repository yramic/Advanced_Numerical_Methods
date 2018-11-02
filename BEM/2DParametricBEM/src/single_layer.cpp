/**
 * \file single_layer.cpp
 * \brief This file declares the functions to evaluate the entries of
 *        Galerkin matrices based on the bilinear form induced by the
 *        Single Layer BIO, using the transformations given in section
 *        1.4.3.4 in the Lecture Notes for Advanced Numerical Methods
 *        for CSE.
 *
 * This File is a part of the 2D-Parametric BEM package
 */

#include "single_layer.hpp"

#include <math.h>
#include <vector>

#include "abstract_bem_space.hpp"
#include "abstract_parametrized_curve.hpp"
#include "discontinuous_space.hpp"
#include "gauleg.hpp"
#include "logweight_quadrature.hpp"
#include "parametrized_mesh.hpp"
#include <Eigen/Dense>

namespace parametricbem2d {
namespace single_layer {
Eigen::MatrixXd InteractionMatrix(const AbstractParametrizedCurve &pi,
                                  const AbstractParametrizedCurve &pi_p,
                                  const AbstractBEMSpace &space,
                                  const unsigned int &N) {
  double tol = 1e-5;

  if (&pi == &pi_p) // Same Panels case
    return ComputeIntegralCoinciding(pi, pi_p, space, N);

  else if (fabs((pi(1) - pi_p(-1)).norm()) < tol ||
           fabs((pi(-1) - pi_p(1)).norm()) < tol) // Adjacent Panels case
    return ComputeIntegralAdjacent(pi, pi_p, space, N);

  else // Disjoint panels case
    return ComputeIntegralGeneral(pi, pi_p, space, N);
}

Eigen::MatrixXd ComputeIntegralCoinciding(const AbstractParametrizedCurve &pi,
                                          const AbstractParametrizedCurve &pi_p,
                                          const AbstractBEMSpace &space,
                                          const unsigned int &N) {
  int Q = space.getQ(); // No. of Reference Shape Functions in trial/test space
  // Interaction matrix with size Q x Q
  Eigen::MatrixXd interaction_matrix(Q, Q);
  // Computing the (i,j)th matrix entry
  for (int i = 0; i < Q; ++i) {
    for (int j = 0; j < Q; ++j) {
      // Lambda expression for functions F and G in eq. 1.4.154
      auto F = [&](double t) { // Function associated with panel pi_p
        return space.evaluateShapeFunction(j, t) * pi_p.Derivative(t).norm();
      };

      auto G = [&](double s) { // Function associated with panel pi
        return space.evaluateShapeFunction(i, s) * pi.Derivative(s).norm();
      };

      // Lambda expression for the 1st integrand in eq. 1.4.159
      auto integrand1 = [&](double s, double t) {
        double tol = 1e-5;
        double s_st;
        if (fabs(s - t) > tol) // Away from singularity
          // Simply evaluating the expression
          s_st = (pi(s) - pi_p(t)).squaredNorm() / (s - t) / (s - t);
        else // Near singularity
          // Using analytic limit for s - > t
          s_st = pi.Derivative(t).squaredNorm();
        return 0.5 * log(s_st) * F(t) * G(s);
      };

      double i1 = 0., i2 = 0.; // The two integrals in eq. 1.4.159

      // Getting Gauss Legendre quadrature weights and points
      Eigen::RowVectorXd weights, points;
      std::tie(points, weights) = gauleg(-1, 1, N);

      // Tensor product quadrature for double 1st integral in eq. 1.4.159
      for (unsigned int i = 0; i < N; ++i) {
        for (unsigned int j = 0; j < N; ++j) {
          i1 += weights(i) * weights(j) * integrand1(points(i), points(j));
        }
      }

      // Lambda expression for 2nd integrand (1.4.159) in transformed
      // coordinates (1.4.163)
      auto integrand2 = [&](double w, double z) {
        return F(0.5 * (w - z)) * G(0.5 * (w + z)) +
               F(0.5 * (w + z)) * G(0.5 * (w - z));
      };

      // Getting log weighted quadrature nodes and weights
      QuadRule logweightQR = getLogWeightQR(2., N);

      // Double loop for 2nd double integral (1.4.163)
      for (unsigned int i = 0; i < N; ++i) {
        // Outer integral evaluated with Log weighted quadrature
        double z = logweightQR.x(i);
        double inner = 0.;
        // Evaluating the inner integral for fixed z
        for (unsigned int j = 0; j < N; ++j) {
          // Scaling Gauss Legendre quadrature nodes to the integral limits
          double w = points(j) * (2 - z);
          inner += weights(j) * integrand2(w, z);
        }
        // Multiplying the integral with appropriate constants for
        // transformation to w from Gauss Legendre nodes
        inner *= (2 - z);
        i2 += logweightQR.w(i) * inner;
      }
      // Filling the matrix entry
      interaction_matrix(i, j) = -1. / (2. * M_PI) * (i1 + 0.5 * i2);
    }
  }
  return interaction_matrix;
}

Eigen::MatrixXd ComputeIntegralAdjacent(const AbstractParametrizedCurve &pi,
                                        const AbstractParametrizedCurve &pi_p,
                                        const AbstractBEMSpace &space,
                                        const unsigned int &N) {
  int Q = space.getQ(); // No. of Reference Shape Functions in trial/test space
  // Interaction matrix with size Q x Q
  Eigen::MatrixXd interaction_matrix(Q, Q);
  // Computing the (i,j)th matrix entry
  for (int i = 0; i < Q; ++i) {
    for (int j = 0; j < Q; ++j) {
      // Panel lengths for local arclength parametrization. Actual values are
      // not required so a length of 1 is used for both the panels
      double length_pi = 1.;   // Length for panel pi
      double length_pi_p = 1.; // Length for panel pi_p

      // when transforming the parametrizations from [-1,1]->\Pi to local
      // arclength parametrizations [0,|\Pi|] -> \Pi, swap is used to ensure
      // that the common point between the panels corresponds to the parameter 0
      // in both arclength parametrizations
      bool swap = (pi(1) != pi_p(-1));

      // Lambda expressions for the functions F,G and D(r,phi) in eq. 1.4.172
      auto F = [&](double t_pr) { // Function associated with panel pi_p
        // Transforming the local arclength parameter to standard parameter
        // range [-1,1] using swap
        double t =
            swap ? 1 - 2 * t_pr / length_pi_p : 2 * t_pr / length_pi_p - 1;
        return space.evaluateShapeFunction(j, t) * pi_p.Derivative(t).norm();
      };

      auto G = [&](double s_pr) { // Function associated with panel pi
        // Transforming the local arclength parameter to standard parameter
        // range [-1,1] using swap
        double s = swap ? 2 * s_pr / length_pi - 1 : 1 - 2 * s_pr / length_pi;
        return space.evaluateShapeFunction(i, s) * pi.Derivative(s).norm();
      };

      auto D_r_phi = [&](double r, double phi) { // eq. 1.4.172
        // Transforming to local arclength parameter range
        double s_pr = r * cos(phi);
        // Transforming to standard parameter range [-1,1] using swap
        double s = swap ? 2 * s_pr / length_pi - 1 : 1 - 2 * s_pr / length_pi;
        // Transforming to local arclength parameter range
        double t_pr = r * sin(phi);
        // Transforming to standard parameter range [-1,1] using swap
        double t =
            swap ? 1 - 2 * t_pr / length_pi_p : 2 * t_pr / length_pi_p - 1;
        if (r != 0) // Away from singularity, simply use the formula
          return (pi(s) - pi_p(t)).squaredNorm() / r / r;
        else // Near singularity, use analytically evaluated limit for r -> 0
          return 1 - sin(2 * phi) * pi.Derivative(s).dot(pi_p.Derivative(t));
      };

      // Getting Gauss Legendre Quadrature weights and nodes
      Eigen::RowVectorXd weights, points;
      std::tie(points, weights) = gauleg(-1, 1, N);

      // The two integrals in eq. 1.4.172 have to be further split into two
      // parts part 1 is where phi goes from 0 to alpha part 2 is where phi goes
      // from alpha to pi/2
      double alpha = atan(length_pi_p / length_pi); // the split point

      // i_IJ -> Integral I, part J
      double i11 = 0., i21 = 0., i12 = 0., i22 = 0.;
      // part 1 (phi from 0 to alpha)
      for (unsigned int i = 0; i < N; ++i) {
        // Transforming gauss quadrature node into phi
        double phi = alpha / 2 * (1 + points(i));
        // Computing inner integral with fixed phi
        // Inner integral for double integral 1, evaluated with Gauss Legendre
        // quadrature
        double inner1 = 0.;
        // Inner integral for double integral 2, evaluated with Log weighted
        // Gauss quadrature
        double inner2 = 0.;
        // Upper limit for inner 'r' integral
        double rmax = length_pi / cos(phi);
        // Evaluating the inner 'r' integral
        for (unsigned int j = 0; j < N; ++j) {
          // Getting Quadrature weights and nodes for Log weighted Gauss
          // quadrature
          QuadRule logweightQR = getLogWeightQR(rmax, N);
          // Evaluating inner2 using Log weighted Gauss quadrature
          double r = logweightQR.x(j);
          inner2 += logweightQR.w(j) * r * F(r * sin(phi)) * G(r * cos(phi));

          // Evaluating inner1 using Gauss Legendre quadrature
          r = rmax / 2 * (1 + points(j));
          inner1 += weights(j) * r * log(D_r_phi(r, phi)) * F(r * sin(phi)) *
                    G(r * cos(phi));
        }
        // Multiplying the integral with appropriate constants for
        // transformation to r from Gauss Legendre nodes
        inner1 *= rmax / 2;
        // Multiplying the integrals with appropriate constants for
        // transformation to phi from Gauss Legendre nodes
        i11 += weights(i) * inner1 * alpha / 2;
        i21 += weights(i) * inner2 * alpha / 2;
      }

      // part 2 (phi from alpha to pi/2)
      for (unsigned int i = 0; i < N; ++i) {
        // Transforming gauss quadrature node into phi (alpha to pi/2)
        double phi =
            points(i) * (M_PI / 2. - alpha) / 2. + (M_PI / 2. + alpha) / 2.;
        // Computing inner integral with fixed phi
        // Inner integral for double integral 1, evaluated with Gauss Legendre
        // quadrature
        double inner1 = 0.;
        // Inner integral for double integral 2, evaluated with Log weighted
        // Gauss quadrature
        double inner2 = 0.;
        // Upper limit for inner 'r' integral
        double rmax = length_pi_p / sin(phi);
        // Evaluating the inner 'r' integral
        for (unsigned int j = 0; j < N; ++j) {
          // Getting Quadrature weights and nodes for Log weighted Gauss
          // quadrature
          QuadRule logweightQR = getLogWeightQR(rmax, N);
          // Evaluating inner2 using Log weighted Gauss quadrature
          double r = logweightQR.x(j);
          inner2 += logweightQR.w(j) * r * F(r * sin(phi)) * G(r * cos(phi));

          // Evaluating inner1 using Gauss Legendre quadrature
          r = rmax / 2 * (1 + points(j));
          inner1 += weights(j) * r * log(D_r_phi(r, phi)) * F(r * sin(phi)) *
                    G(r * cos(phi));
        }
        // Multiplying the integral with appropriate constants for
        // transformation to r from Gauss Legendre quadrature nodes
        inner1 *= rmax / 2;
        // Multiplying the integrals with appropriate constants for
        // transformation to phi from Gauss Legendre quadrature nodes
        i12 += weights(i) * inner1 * (M_PI / 2. - alpha) / 2.;
        i22 += weights(i) * inner2 * (M_PI / 2. - alpha) / 2.;
      }
      // Summing up the parts to get the final integral
      double integral = 0.5 * (i11 + i12) + (i21 + i22);
      // Multiplying the integral with appropriate constants for transformation
      // to local arclength variables
      integral *= 4 / length_pi / length_pi_p;
      // Filling up the matrix entry
      interaction_matrix(i, j) = -1 / (2 * M_PI) * integral;
    }
  }
  return interaction_matrix;
}

Eigen::MatrixXd ComputeIntegralGeneral(const AbstractParametrizedCurve &pi,
                                       const AbstractParametrizedCurve &pi_p,
                                       const AbstractBEMSpace &space,
                                       const unsigned int &N) {
  int Q = space.getQ(); // No. of Reference Shape Functions in trial/test space
  // Interaction matrix with size Q x Q
  Eigen::MatrixXd interaction_matrix(Q, Q);
  // Computing the (i,j)th matrix entry
  for (int i = 0; i < Q; ++i) {
    for (int j = 0; j < Q; ++j) {
      // Lambda expression for functions F and G in eq. 1.4.154
      auto F = [&](double t) { // Function associated with panel pi_p
        return space.evaluateShapeFunction(j, t) * pi_p.Derivative(t).norm();
      };

      auto G = [&](double s) { // Function associated with panel pi
        return space.evaluateShapeFunction(i, s) * pi.Derivative(s).norm();
      };

      double integral = 0.;

      // Getting Gauss Legendre quadrature weights and points
      Eigen::RowVectorXd weights, points;
      std::tie(points, weights) = gauleg(-1, 1, N);

      // Tensor product quadrature rule
      for (unsigned int i = 0; i < N; ++i) {
        for (unsigned int j = 0; j < N; ++j) {
          double s = points(i);
          double t = points(j);
          integral += weights(i) * weights(j) * log((pi(s) - pi_p(t)).norm()) *
                      F(t) * G(s);
        }
      }
      // Filling up the matrix entry
      interaction_matrix(i, j) = -1 / (2 * M_PI) * integral;
    }
  }
  return interaction_matrix;
}

Eigen::MatrixXd SingleLayerMatrix(const ParametrizedMesh mesh,
                                  const AbstractBEMSpace &space,
                                  const unsigned int &N) {
  // Getting the number of panels in the mesh
  unsigned int numpanels = mesh.getNumPanels();
  // Getting dimensions of trial/test space
  unsigned int dims = space.getSpaceDim(numpanels);
  // Getting the panels from the mesh
  PanelVector panels = mesh.getPanels();
  // Getting the number of local shape functions in the trial/test space
  unsigned int Q = space.getQ();
  // Initializing the Galerkin matrix with zeros
  Eigen::MatrixXd output = Eigen::MatrixXd::Zero(dims, dims);
  // Panel oriented assembly
  for (unsigned int i = 0; i < numpanels; ++i) {
    for (unsigned int j = 0; j < numpanels; ++j) {
      // Getting the interaction matrix for the pair of panels i and j
      Eigen::MatrixXd interaction_matrix =
          InteractionMatrix(*panels[i], *panels[j], space, N);
      // Local to global mapping of the elements in interaction matrix
      for (unsigned int I = 0; I < Q; ++I) {
        for (unsigned int J = 0; J < Q; ++J) {
          int II = space.LocGlobMap(I + 1, i + 1, numpanels) - 1;
          int JJ = space.LocGlobMap(J + 1, j + 1, numpanels) - 1;
          // Filling the Galerkin matrix entries
          output(II, JJ) += interaction_matrix(I, J);
        }
      }
    }
  }
  return output;
}

} // namespace single_layer
} // namespace parametricbem2d
