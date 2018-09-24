/**
 * \file abstract_bem_space.hpp
 * \brief This file declares an abstract class representing a BEM Space
 *
 * This File is a part of the 2D-Parametric BEM package
 */

#ifndef CONTSPACEHPP
#define CONTSPACEHPP

#include "abstract_bem_space.hpp"

#include <vector>
#include <utility>

namespace parametricbem2d {
  /**
   * \class AbstractBEMSpace
   * \brief This abstract class declares the interface for a BEM Space
   */
  template <unsigned int N>
  class ContinuousSpace : public AbstractBEMSpace {
  public:
    ContinuousSpace() {
      std::cout << "Error! No specialization defined" << std::endl;
    }
  }; // class AbstractBEMSpace

  // Specialization of the class for N = 0
  template <>
  class ContinuousSpace<0> : public AbstractBEMSpace{
  public:
    static double b1(double t) {
      assert(IsWithinParameterRange(t));
      return 1.;
    }
    // Constructor
    ContinuousSpace() {
      q = 1;
      referenceshapefunctions.push_back(b1);
    }
  };

  // Specialization of the class for N = 1
  template <>
  class ContinuousSpace<1> : public AbstractBEMSpace{
  public:
    static double b1(double t) {
      assert(IsWithinParameterRange(t));
      return 0.5*(t+1);
    }
    static double b2(double t) {
      assert(IsWithinParameterRange(t));
      return 0.5*(1-t);
    }
    // Constructor
    ContinuousSpace() {
      q = 2;
      referenceshapefunctions.push_back(b1);
      referenceshapefunctions.push_back(b2);
    }
  };

  // Specialization of the class for N = 1
  template <>
  class ContinuousSpace<2> : public AbstractBEMSpace{
  public:
    static double b1(double t) {
      assert(IsWithinParameterRange(t));
      return 0.5*(t+1);
    }
    static double b2(double t) {
      assert(IsWithinParameterRange(t));
      return 0.5*(1-t);
    }
    static double b3(double t) {
      assert(IsWithinParameterRange(t));
      return 1-t*t;
    }
    // Constructor
    ContinuousSpace() {
      q = 3;
      referenceshapefunctions.push_back(b1);
      referenceshapefunctions.push_back(b2);
      referenceshapefunctions.push_back(b3);
    }
  };
} // namespace parametricbem2d

#endif//CONTSPACEHPP
