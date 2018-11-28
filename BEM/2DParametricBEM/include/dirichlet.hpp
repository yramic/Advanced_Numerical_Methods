/**
 * \file dirichlet.hpp
 * \brief This file declares the functions to evaluate the entries of
 *        Galerkin matrices based on the bilinear form induced by the
 *        Hypersingular BIO, using the transformations given in section
 *        1.4.3.4 in the Lecture Notes for Advanced Numerical Methods
 *        for CSE.
 *
 * This File is a part of the 2D-Parametric BEM package
 */

#ifndef DIRICHLETHPP
#define DIRICHLETHPP

#include <Eigen/Dense>

#include "abstract_bem_space.hpp"
#include "abstract_parametrized_curve.hpp"
#include "parametrized_mesh.hpp"

namespace parametricbem2d {
/**
 * This namespace contains all the functions for evaluating the Single Layer
 * Galerkin Matrix using quadrature and panel oriented assembly
 */
 inline Eigen::MatrixXd MassMatrix(const ParametrizedMesh &mesh, const AbstractBEMSpace& x_space,const AbstractBEMSpace& y_space, unsigned order) {
   unsigned qx = x_space.getQ();
   unsigned qy = y_space.getQ();
   PanelVector panels = mesh.getPanels();
   unsigned numpanels = mesh.getNumPanels();
   unsigned int rows = x_space.getSpaceDim(numpanels);
   unsigned int cols = y_space.getSpaceDim(numpanels);
   Eigen::MatrixXd output = Eigen::MatrixXd::Zero(rows, cols);
   Eigen::MatrixXd interaction(qx,qy);
   for (unsigned panel = 0 ; panel < numpanels ; ++panel) {
     for (unsigned i = 0 ; i < qx ; ++i) {
       for (unsigned j = 0 ; j < qy ; ++j) {
         std::function<double(double)> integrand = [&] (double x) {
           return x_space.evaluateShapeFunction(i,x) * y_space.evaluateShapeFunction(j,x) * panels[panel]->Derivative(x).norm();
         };
         interaction(i,j) = ComputeIntegral(integrand,-1,1,order);
       }
     }
     // Local to global mapping of the elements in interaction
     for (unsigned int I = 0; I < qx; ++I) {
       for (unsigned int J = 0; J < qy; ++J) {
         int II = x_space.LocGlobMap(I + 1, panel + 1, numpanels) - 1;
         int JJ = y_space.LocGlobMap(J + 1, panel + 1, numpanels) - 1;
         // Filling the Galerkin matrix entries
         output(II, JJ) += interaction(I, J);
       }
     }
   }
   return output;
 }
/*namespace dirichlet_bvp {
  namespace direct_first_kind {
    Eigen::MatrixXd solve(const ParametrizedMesh &mesh, std::function<double(double)> g) {
      DiscontinuousSpace<0> trial_space; // Same as test space
      DiscontinuousSpace<0> test_space;
      ContinuousSpace<1> g_interpol_space;
      Eigen::MatrixXd V = single_layer::GalerkinMatrix(mesh,trial_space,32);
      Eigen::MatrixXd K = double_layer::GalerkinMatrix(mesh,g_interpol_space,test_space,32);
      Eigen::MatrixXd M = MassMatrix(mesh,trial_space,test_space,32);
    }
  }

  namespace direct_second_kind {

  }

  namespace direct_second_kind {

  }

  namespace direct_second_kind {

  }*/

} // namespace parametricbem2d

#endif // DIRICHLETHPP
