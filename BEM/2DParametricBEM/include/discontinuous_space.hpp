/**
 * \file discontinuous_space.hpp
 * \brief This file declares a templated class inherited from AbstractBEMSpace
 *        and represents discontinuous spaces of the form \f$S^{-1}_{p}\f$ (eq.1.4.22)
 *        using full template specialization.
 *
 * This File is a part of the 2D-Parametric BEM package
 */

#ifndef DISCONTSPACEHPP
#define DISCONTSPACEHPP

#include "abstract_bem_space.hpp"

#include <iostream>
#include <vector>
#include <utility>

namespace parametricbem2d {
  /**
   * \class DiscontinuousSpace
   * \brief This template class inherits from the class AbstractBEMSpace
   *        . This class implements the BEM spaces of the form \f$S^{-1}_{p}\f$
   *        defined in eq. .The class is implemented through full template
   *        specialization for different values of p.
   */
  template <unsigned int p>
  class DiscontinuousSpace : public AbstractBEMSpace {
  public:
    DiscontinuousSpace() {
      std::cout << "Error! No specialization defined for p = " << p << std::endl;
    }
  };

  /**
   * \brief This is a specialization of the templated class for p = 0.
   *        This class represents the space \f$S^{-1}_{0}\f$
   */
  template <>
  class DiscontinuousSpace<0> : public AbstractBEMSpace{
  public:
    // Basis function
    static double b1(double t) {
      assert(IsWithinParameterRange(t));
      return 1.;
    }
    // Local to Global Map
    static unsigned int map(unsigned int q, unsigned int n, unsigned int N) {
      assert(q<=1 && n<=N);
      return n;
    }
    // Space Dimensions
    unsigned int getSpaceDim(unsigned int numpanels) const{
      return numpanels*q_;
    }
    // Constructor
    DiscontinuousSpace() {
      q_ = 1;
      referenceshapefunctions_.push_back(b1);
      locglobmap_ = map;
    }
  };

  /**
   * \brief This is a specialization of the templated class for p = 1.
   *        This class represents the space \f$S^{-1}_{1}\f$
   */
  template <>
  class DiscontinuousSpace<1> : public AbstractBEMSpace{
  public:
    // Basis function 1
    static double b1(double t) {
      assert(IsWithinParameterRange(t));
      return 0.5;
    }
    // Basis function 2
    static double b2(double t) {
      assert(IsWithinParameterRange(t));
      return 0.5*t;
    }
    // Local to Global Map
    static unsigned int map(unsigned int q, unsigned int n, unsigned int N) {
      assert(q<=2 && n<=N);
      if (q==1)
        return n;
      else
        return N+n;
    }
    // Space Dimensions
    unsigned int getSpaceDim(unsigned int numpanels) const{
      return numpanels*q_;
    }
    // Constructor
    DiscontinuousSpace() {
      q_ = 2;
      referenceshapefunctions_.push_back(b1);
      referenceshapefunctions_.push_back(b2);
      locglobmap_ = map;
    }
  };
} // namespace parametricbem2d

#endif//DISCONTSPACEHPP
