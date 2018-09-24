/**
 * \file abstract_bem_space.hpp
 * \brief This file declares an abstract class representing a BEM Space
 *
 * This File is a part of the 2D-Parametric BEM package
 */

#ifndef DISCONTSPACEHPP
#define DISCONTSPACEHPP

#include "abstract_bem_space.hpp"

#include <vector>
#include <utility>

namespace parametricbem2d {
  /**
   * \class AbstractBEMSpace
   * \brief This abstract class declares the interface for a BEM Space
   */
  template <unsigned int N>
  class DiscontinuousSpace : public AbstractBEMSpace {
  public:
    DiscontinuousSpace():p(N) {
      std::cout << "Error! No specialization defined" << std::endl;
    }
  private:
    int p;
  }; // class AbstractBEMSpace

  // Specialization of the class for N = 0
  template <>
  class DiscontinuousSpace<0> : public AbstractBEMSpace{
  public:
    static double b1(double t) {
      assert(IsWithinParameterRange(t));
      return 1.;
    }
    // Constructor
    DiscontinuousSpace() {
      q = 1;
      referenceshapefunctions.push_back(b1);
    }
  };

  // Specialization of the class for N = 1
  template <>
  class DiscontinuousSpace<1> : public AbstractBEMSpace{
  public:
    static double b1(double t) {
      assert(IsWithinParameterRange(t));
      return 0.5;
    }
    static double b2(double t) {
      assert(IsWithinParameterRange(t));
      return 0.5*t;
    }
    // Constructor
    DiscontinuousSpace() {
      q = 2;
      referenceshapefunctions.push_back(b1);
      referenceshapefunctions.push_back(b2);
    }
  };
} // namespace parametricbem2d

#endif//DISCONTSPACEHPP
