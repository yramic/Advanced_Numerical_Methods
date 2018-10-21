/**
 * \file double_layer.cpp
 * \brief This file declares the functions to evaluate the entries of
 *        Galerkin matrices based on the bilinear form induced by the
 *        Double Layer BIO, using the transformations given in section
 *        1.4.3.4 in the Lecture Notes for Advanced Numerical Methods
 *        for CSE.
 *
 * This File is a part of the 2D-Parametric BEM package
 */

#include "double_layer.hpp"

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
namespace double_layer {
Eigen::MatrixXd InteractionMatrix(const AbstractParametrizedCurve &pi,
                                  const AbstractParametrizedCurve &pi_p,
                                  const AbstractBEMSpace &trial_space,
                                  const AbstractBEMSpace &test_space,
                                  const unsigned int &N) {
  double tol = 1e-5;
  if (&pi == &pi_p) // Same Panels
    return ComputeIntegralCoinciding(pi, pi_p, trial_space, test_space, N);
  else if (fabs((pi(1) - pi_p(-1)).norm()) < tol ||
           fabs((pi(-1) - pi_p(1)).norm()) < tol) // Adjacent Panels
    return ComputeIntegralAdjacent(pi, pi_p, trial_space, test_space, N);
  else
    return ComputeIntegralGeneral(pi, pi_p, trial_space, test_space, N);
}

Eigen::MatrixXd ComputeIntegralCoinciding(const AbstractParametrizedCurve &pi,
                                          const AbstractParametrizedCurve &pi_p,
                                          const AbstractBEMSpace &trial_space,
                                          const AbstractBEMSpace &test_space,
                                          const unsigned int &N) {
  //
  int Qtrial =
      trial_space.getQ(); // The number of Reference Shape Functions in space
  int Qtest =
      test_space.getQ(); // The number of Reference Shape Functions in space
  Eigen::MatrixXd interaction_matrix(Qtest, Qtrial);
  for (int i = 0; i < Qtest; ++i) {
    for (int j = 0; j < Qtrial; ++j) {
      // Lambda expression for functions F and G
      auto F = [&](double t) {
        return trial_space.evaluateShapeFunction(j, t) *
               pi_p.Derivative(t).norm();
      };

      auto G = [&](double s) {
        return test_space.evaluateShapeFunction(i, s) * pi.Derivative(s).norm();
      };

      // Lambda expression for the integrand in eq. 1.4.174
      auto integrand = [&](double s, double t) {
        double tol = 1e-5;
        double k;
        Eigen::Vector2d tangent = pi_p.Derivative(t);
        // Normal vector assuming that the curve is counter clockwise
        Eigen::Vector2d normal;
        normal << -tangent(1), tangent(0);
        normal = normal / normal.norm();
        if (fabs(s - t) > tol) {
          k = (pi(s) - pi_p(t)).dot(normal) / (pi(s) - pi_p(t)).squaredNorm();
          // std::cout<< "s,t far, dot = " << k << std::endl;
        }

        else {
          // Limit evaluated analytically for s -> t
          k = 0.5 * pi.DoubleDerivative(t).dot(normal) /
              pi.Derivative(t).squaredNorm();
          // std::cout<< "s,t near, dot = " << k << std::endl;
        }

        return k * F(t) * G(s);
      };

      double integral = 0;

      // Getting quadrature weights and points
      Eigen::RowVectorXd weights, points;
      std::tie(points, weights) = gauleg(-1, 1, N);

      // Double sum for double integral in eq. 1.4.174
      for (unsigned int i = 0; i < N; ++i) {
        for (unsigned int j = 0; j < N; ++j) {
          integral += weights(i) * weights(j) * integrand(points(i), points(j));
        }
      }
      interaction_matrix(i, j) = -1. / (2. * M_PI) * integral;
    }
  }
  return interaction_matrix;
}

Eigen::MatrixXd ComputeIntegralAdjacent(const AbstractParametrizedCurve &pi,
                                        const AbstractParametrizedCurve &pi_p,
                                        const AbstractBEMSpace &trial_space,
                                        const AbstractBEMSpace &test_space,
                                        const unsigned int &N) {
  //
  int Qtrial =
      trial_space.getQ(); // The number of Reference Shape Functions in space
  int Qtest =
      test_space.getQ(); // The number of Reference Shape Functions in space
  Eigen::MatrixXd interaction_matrix(Qtest, Qtrial);
  for (int i = 0; i < Qtest; ++i) {
    for (int j = 0; j < Qtrial; ++j) {
      // Getting the panel lengths for local arclength parametrization
      double length_pi = 1.; // Panel length for local arclength parametrization
      double length_pi_p =
          1.; // Panel length for local arclength parametrization

      // swap ensures when transforming the parametrizations from [-1,1]->\Pi to
      // local arclength parametrizations [0,|\Pi|] -> \Pi, the common point
      // between the panels corresponds to the parameter 0 in arclength
      // parametrizations
      bool swap = (pi(1) != pi_p(-1));

      // Lambda expressions for the functions F,G and D(r,phi)
      auto F = [&](double t_pr) { // Function associated with panel pi_p
        double t =
            swap ? 1 - 2 * t_pr / length_pi_p : 2 * t_pr / length_pi_p - 1;
        return trial_space.evaluateShapeFunction(j, t) *
               pi_p.Derivative(t).norm();
      };

      auto G = [&](double s_pr) { // Function associated with panel pi
        double s = swap ? 2 * s_pr / length_pi - 1 : 1 - 2 * s_pr / length_pi;
        return test_space.evaluateShapeFunction(i, s) * pi.Derivative(s).norm();
      };

      auto integrand = [&](double r, double phi) { // eq. 1.4.172
        double s_pr = r * cos(phi);
        double s = swap ? 2 * s_pr / length_pi - 1 : 1 - 2 * s_pr / length_pi;
        double t_pr = r * sin(phi);
        double t =
            swap ? 1 - 2 * t_pr / length_pi_p : 2 * t_pr / length_pi_p - 1;
        Eigen::Vector2d tangent = pi_p.Derivative(t);
        // Normal vector assuming that the curve is counter clockwise
        Eigen::Vector2d normal;
        normal << -tangent(1), tangent(0);
        normal = normal / normal.norm();
        if (r > 1e-5)
          return r * (pi(s) - pi_p(t)).dot(normal) /
                 (pi(s) - pi_p(t)).squaredNorm();
        else {
          double s0 = swap ? -1 : 1;
          double t0 = swap ? 1 : -1;
          Eigen::Vector2d b_r_phi =
              pi.Derivative(s0) * cos(phi) - pi_p.Derivative(t0) * sin(phi);
          return b_r_phi.dot(normal) / b_r_phi.squaredNorm();
        }
      };

      // Getting Gauss Quadrature weights and nodes
      Eigen::RowVectorXd weights, points;
      std::tie(points, weights) = gauleg(-1, 1, N);

      // The integral is split into two parts according to eq. 1.4.178
      // part 1 is where phi goes from 0 to alpha
      // part 2 is where phi goes from alpha to pi/2
      double alpha = atan(length_pi_p / length_pi);

      // i_J -> part J
      double i1 = 0., i2 = 0.;
      // part 1
      for (unsigned int i = 0; i < N; ++i) {
        // Transforming gauss quadrature abcissa into phi (0 to alpha)
        double phi = alpha / 2 * (1 + points(i));
        // Computing inner integral with fixed phi
        double inner = 0.; // inner integral
        double rmax =
            length_pi / cos(phi); // Upper limit for inner 'r' integral
        for (unsigned int j = 0; j < N; ++j) {
          double r =
              rmax / 2 * (1 + points(j)); // standard gauss quadrature for inner
          inner += weights(j) * integrand(r, phi) * F(r * sin(phi)) *
                   G(r * cos(phi));
        }
        inner *= rmax / 2;
        i1 += weights(i) * inner * alpha / 2;
      }

      // part 2 (phi from alpha to pi/2)
      for (unsigned int i = 0; i < N; ++i) {
        // Transforming gauss quadrature abcissa into phi (alpha to pi/2)
        double phi =
            points(i) * (M_PI / 2. - alpha) / 2. + (M_PI / 2. + alpha) / 2.;
        double inner = 0.; // inner integral
        double rmax =
            length_pi_p / sin(phi); // Upper limit for inner 'r' integral
        for (unsigned int j = 0; j < N; ++j) {
          double r = rmax / 2 *
                     (1 + points(j)); // for inner1 standard gauss quadrature
          inner += weights(j) * integrand(r, phi) * F(r * sin(phi)) *
                   G(r * cos(phi));
        }
        inner *= rmax / 2;
        i2 += weights(i) * inner * (M_PI / 2. - alpha) / 2.;
      }
      double integral = i1 + i2;
      integral *= 4 / length_pi / length_pi_p;
      interaction_matrix(i, j) = -1 / (2 * M_PI) * integral;
    }
  }
  return interaction_matrix;
}

Eigen::MatrixXd ComputeIntegralGeneral(const AbstractParametrizedCurve &pi,
                                       const AbstractParametrizedCurve &pi_p,
                                       const AbstractBEMSpace &trial_space,
                                       const AbstractBEMSpace &test_space,
                                       const unsigned int &N) {
  //
  int Qtrial =
      trial_space.getQ(); // The number of Reference Shape Functions in space
  int Qtest =
      test_space.getQ(); // The number of Reference Shape Functions in space
  Eigen::MatrixXd interaction_matrix(Qtest, Qtrial);
  for (int i = 0; i < Qtest; ++i) {
    for (int j = 0; j < Qtrial; ++j) {
      // Lambda expression for functions F and G
      auto F = [&](double t) {
        return trial_space.evaluateShapeFunction(j, t) *
               pi_p.Derivative(t).norm();
      };

      auto G = [&](double s) {
        return test_space.evaluateShapeFunction(i, s) * pi.Derivative(s).norm();
      };

      auto integrand = [&](double s, double t) { // eq. 1.4.172
        Eigen::Vector2d tangent = pi_p.Derivative(t);
        // Normal vector assuming that the curve is counter clockwise
        Eigen::Vector2d normal;
        normal << -tangent(1), tangent(0);
        normal = normal / normal.norm();
        return (pi(s) - pi_p(t)).dot(normal) / (pi(s) - pi_p(t)).squaredNorm();
      };

      double integral = 0.;

      // Getting quadrature weights and points
      Eigen::RowVectorXd weights, points;
      std::tie(points, weights) = gauleg(-1, 1, N);

      // Double sum for double integral
      for (unsigned int i = 0; i < N; ++i) {
        for (unsigned int j = 0; j < N; ++j) {
          double s = points(i);
          double t = points(j);
          integral += weights(i) * weights(j) * integrand(s, t) * F(t) * G(s);
        }
      }
      interaction_matrix(i, j) = -1 / (2 * M_PI) * integral;
    }
  }
  return interaction_matrix;
}

Eigen::MatrixXd DoubleLayerMatrix(const ParametrizedMesh mesh,
                                  const AbstractBEMSpace &trial_space,
                                  const AbstractBEMSpace &test_space,
                                  const unsigned int &N) {
  unsigned int numpanels = mesh.getNumPanels();
  unsigned int rows = test_space.getSpaceDim(numpanels);
  unsigned int cols = trial_space.getSpaceDim(numpanels);
  PanelVector panels = mesh.getPanels();
  unsigned int Qtest = test_space.getQ();
  unsigned int Qtrial = trial_space.getQ();
  Eigen::MatrixXd output = Eigen::MatrixXd::Zero(rows, cols);

  for (unsigned int i = 0; i < numpanels; ++i) {
    for (unsigned int j = 0; j < numpanels; ++j) {
      Eigen::MatrixXd interaction_matrix =
          InteractionMatrix(*panels[i], *panels[j], trial_space, test_space, N);
      // Local to global mapping
      for (unsigned int I = 0; I < Qtest; ++I) {
        for (unsigned int J = 0; J < Qtrial; ++J) {
          int II = test_space.LocGlobMap(I + 1, i + 1, numpanels) - 1;
          int JJ = trial_space.LocGlobMap(J + 1, j + 1, numpanels) - 1;
          output(II, JJ) += interaction_matrix(I, J);
        }
      }
    }
  }
  return output;
}

} // namespace double_layer
} // namespace parametricbem2d
