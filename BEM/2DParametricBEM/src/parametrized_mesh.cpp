/**
 * \file parametrized_mesh.cpp
 * \brief This file defines the class for representing a mesh comprised
 *        of panels represented by parametrized curves.
 *
 *  This File is a part of the 2D-Parametric BEM package
 */

#include "parametrized_mesh.hpp"

#include <utility>

#include <Eigen/Dense>

namespace parametricbem2d {

  ParametrizedMesh::ParametrizedMesh(PanelVector panels):panels_(panels) {
   unsigned int N = getNumPanels();
   // Checking consistency of the mesh
   for (unsigned int i = 0 ; i < N ; ++i)
     assert(fabs(panels_[i%N]->(1)-panels_[(i+1)%N]->(-1)) < 1e-5);
  }

  PanelVector ParametrizedMesh::getPanels() const{
   return panels_;
  }

  unsigned int ParametrizedMesh::getNumPanels() const {
   return panels_.size();
  }

  Eigen::Vector2d ParametrizedMesh::getVertex(unsigned int i) const {
   assert(i<getNumPanels());
   return panels_[i]->(-1);
  }

} // namespace parametricbem2d
