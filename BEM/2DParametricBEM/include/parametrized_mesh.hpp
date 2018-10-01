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
   * \brief This abstract class provides the interface for a parametric
   *        curve. The interface follows the one outlined in
   *        Section 1.4.3.4 in the lecture material for
   *        Advanced Numerical Methods for CSE.
   */

  //using PanelVector = std::vector<AbstractParametrizedCurve>;

  class ParametrizedMesh {
    public:
      using PanelIterator = PanelVector::iterator;

     /**
      * This function is used for querying the parameter interval.
      * The standard parameter interval [-1,1] is used and it can't
      * be overriden in the inherited classes as the function is non
      * virtual. The function is declared static as it is independent
      * of the concrete object.
      *
      * @return A std::pair<double,double> object containing
      *          the valid parameter range for parametrization
      *          that is [-1,1]
      */
     ParametrizedMesh(PanelVector panels);

     /**
      * This function is used for retrieving the derivative \f$\dot{\gamma}\f$(t)
      * at a point on the parametrized curve \f$\gamma\f$(t).
      * This is a pure virtual function which has to be implemented
      * in the inherited classes.
      *
      * @param t Parameter value of Double type
      *          at which the derivative for
      *          parametrization is evaluated
      * @return A 2-D vector of type Eigen::Vector2d
      *         containing the parametrized
      *         derivative at point 't'
      */
     PanelVector getPanels() const;

     /**
      * This function is used for checking whether a value t is within the
      * valid parameter range. This function is non virtual to prevent it
      * from being overriden as the parameter interval is fixed. It is
      * declared static because the check is independent of the concrete
      * implementation.
      *
      * @param t The value to be checked
      * @return boolean indicating result of the performed check
      */
     unsigned int getNumPanels() const;

     Eigen::Vector2d getVertex(unsigned int i) const;

   private:
     const PanelVector panels_;
  }; // class ParametrizedMesh
} // namespace parametricbem2d

#endif //PARAMETRIZEDMESHHPP
