/**
 * \file neumann.hpp
 * \brief This file defines various solvers to solve a Neumann Boundary Value
 *        problem of the form given in \f$\eqref{eq:neubvp}\f$
 *
 * This File is a part of the 2D-Parametric BEM package
 */

#ifndef NEUMANNHPP
#define NEUMANNHPP

#include <Eigen/Dense>

#include "abstract_bem_space.hpp"
#include "abstract_parametrized_curve.hpp"
#include "parametrized_mesh.hpp"

namespace parametricbem2d {
/**
 * This function is used to evaluate vectors appearing in augmented variational
 * formulations when the vanishing mean condition is not compatible with
 * standard bases. An entry of this vector is given by
 * \f$ C_{i} = \int_{\Gamma} b^{i}(x) dS(x) \f$
 *
 * @param mesh Parametrized mesh representing the boundary \f$\Gamma\f$.
 * @param space Trial space.
 * @param order The order for gauss/log-weighted quadrature
 * @return An Eigen::VectorXd
 */
inline Eigen::VectorXd MassVector(const ParametrizedMesh &mesh,
                                  const AbstractBEMSpace &space,
                                  unsigned order) {
  // Getting the number of reference shape functions in the space
  unsigned q = space.getQ();
  // Getting the panels
  PanelVector panels = mesh.getPanels();
  unsigned numpanels = mesh.getNumPanels();
  // Getting the space dimensions to fix vector sizes
  unsigned int rows = space.getSpaceDim(numpanels);
  Eigen::VectorXd output = Eigen::VectorXd::Zero(rows);
  // Vector to store evaluation in reference coordinates
  Eigen::VectorXd local(q);
  for (unsigned panel = 0; panel < numpanels; ++panel) {
    for (unsigned i = 0; i < q; ++i) {
      // Transformed integrand (in reference coordinates)
      std::function<double(double)> integrand = [&](double x) {
        return space.evaluateShapeFunction(i, x) *
               panels[panel]->Derivative(x).norm();
      };
      // Evaluating using Gauss quadrature
      local(i) = ComputeIntegral(integrand, -1, 1, order);
    }
    // Local to global mapping of the elements
    for (unsigned int I = 0; I < q; ++I) {
      int II = space.LocGlobMap(I + 1, panel + 1, numpanels) - 1;
      // Filling the vector entries
      output(II) += local(I);
    }
  }
  return output;
}

/**
 * This namespace contains all the solvers for Neumann bvp of the form
 * \f$\eqref{eq:neubvp}\f$. For different methods, the outputs mean something
 * else. Check the solver's documentation for further details.
 */
namespace neumann_bvp {
/**
 * This namespace contains the solver using the direct first kind method which
 * has the variational formulation as given in \f$\eqref{eq:aWdir}\f$. The
 * Solver used the lowest order BEM spaces for computation.
 */
namespace direct_first_kind {
/**
 * This function is used to solve the Neumann boundary value problem given
 * in \f$\eqref{eq:neubvp}\f$ using the variational formulation given in
 * \f$\eqref{eq:aWdir}\f$. The function outputs a vector of estimated Dirichlet
 * trace of the solution u.
 *
 * @param mesh Parametrized mesh representing the boundary \f$\Gamma\f$.
 * @param g Neumann boundary condition in 2D using a function of the form
 *          double(double,double)
 * @param order The order for gauss/log-weighted quadrature
 * @return An Eigen::VectorXd type representing the Dirichlet trace of the
 * solution u
 */
Eigen::VectorXd solve(const ParametrizedMesh &mesh,
                      std::function<double(double, double)> Tn,
                      unsigned order) {
  // Same trial and test spaces
  ContinuousSpace<1> trial_space;
  ContinuousSpace<1> test_space;
  // Space used for interpolation of Neumann data
  DiscontinuousSpace<0>
      Tn_interpol_space;
  // Computing W matrix
  Eigen::MatrixXd W = hypersingular::GalerkinMatrix(mesh, trial_space, order);
  // Computing K matrix
  Eigen::MatrixXd K =
      double_layer::GalerkinMatrix(mesh, Tn_interpol_space, test_space, order);
  // Computing mass matrix
  Eigen::MatrixXd M = MassMatrix(mesh, test_space, Tn_interpol_space, order);
  // Midpoint interpolation of Neumann data
  unsigned numpanels = mesh.getNumPanels();
  Eigen::VectorXd Tn_N(numpanels);
  using PanelVector = parametricbem2d::PanelVector;
  PanelVector panels = mesh.getPanels();
  for (unsigned i = 0; i < numpanels; ++i) {
    Eigen::Vector2d pt = panels[i]->operator()(0.);
    Tn_N(i) = Tn(pt(0), pt(1));
  }
  // The vector c used in augmented formulation
  Eigen::VectorXd c(W.rows() + 1);
  c << MassVector(mesh, test_space, order), 0;
  // Building the augmented lhs matrix for solving.
  Eigen::MatrixXd lhs = Eigen::MatrixXd::Zero(K.rows() + 1, K.cols() + 1);
  // Filling in different blocks
  lhs.block(0, 0, K.rows(), K.cols()) = K;
  lhs.row(K.rows()) = c.transpose();
  lhs.col(K.cols()) = c;
  // Building the augmented rhs matrix for solving
  Eigen::MatrixXd rhs = Eigen::MatrixXd::Zero(K.rows() + 1, K.cols());
  lhs.block(0, 0, K.rows(), K.cols()) = (0.5 * M - K).transpose();
  // Build rhs vector for solving
  Eigen::VectorXd rhs_vector = rhs * Tn_N;
  // Solving for coefficients
  Eigen::VectorXd sol = lhs.lu().solve(rhs_vector);
  return sol;
}
} // namespace direct_first_kind

namespace direct_second_kind {} // namespace direct_second_kind

namespace direct_second_kind {}

namespace direct_second_kind {}
} // namespace neumann_bvp
} // namespace parametricbem2d

#endif // NEUMANNHPP
