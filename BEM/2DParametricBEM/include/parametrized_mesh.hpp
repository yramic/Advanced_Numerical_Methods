/**
 * \file parametrized_mesh.hpp
 * \brief This file declares a class for representing a mesh comprised
 *        of panels represented by parametrized curves.
 *
 *  This File is a part of the 2D-Parametric BEM package
 */

#ifndef PARAMETRIZEDMESHHPP
#define PARAMETRIZEDMESHHPP

#include <utility>
#include <iterator>

#include <Eigen/Dense>
#include "abstract_parametrized_curve.hpp"

namespace parametricbem2d {
  /**
   * \class ParametrizedMesh
   * \brief This class is useful for combining parametrized panels into a
   *        parametrized mesh which is further used in assembly of Galerkin
   *        matrices using the parametric BEM approach.
   */

  class ParametrizedMesh {
    public:
     /**
      * Constructor using a PanelVector object which contains the component
      * panels in parametrized form for the parametrized mesh object.
      */
     ParametrizedMesh(PanelVector panels);

     /**
      * This function is used for retrieving a PanelVector containing all the
      * parametrized panels in the parametrized mesh.
      *
      * @return A PanelVector containing all the parametrized panels in the mesh
      */
     PanelVector getPanels() const;

     /**
      * This function is used for getting the number of panels in the
      * parametrized mesh
      *
      * @return number of panels in the mesh
      */
     unsigned int getNumPanels() const;

     /**
      * This function is used for getting the ith vertex in the parametrized mesh
      *
      * @return ith vertex in the mesh as Eigen::Vector2d
      */
     Eigen::Vector2d getVertex(unsigned int i) const;

   private:
     /**
      * Private const field for the PanelVector of the mesh
      */
     const PanelVector panels_;
  }; // class ParametrizedMesh
} // namespace parametricbem2d

#endif //PARAMETRIZEDMESHHPP
