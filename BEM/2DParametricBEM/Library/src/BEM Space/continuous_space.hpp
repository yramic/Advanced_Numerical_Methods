/**
 * \file continuous_space.hpp
 * \brief This file declares a templated class inherited from AbstractBEMSpace
 *        and represents discontinuous spaces of the form \f$S^{0}_{p}\f$ (eq. 1.4.21)
 *        using full template specialization.
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
   * \class ContinuousSpace
   * \brief This template class inherits from the class AbstractBEMSpace
   *        . This class implements the BEM spaces of the form \f$S^{0}_{p}\f$
   *        defined in eq. .The class is implemented through full template
   *        specialization for different values of p.
   */
  template <unsigned int p>
  class ContinuousSpace : public AbstractBEMSpace {
  public:
    ContinuousSpace() {
      std::cout << "Error! No specialization defined for p = " << p << std::endl;
    }
  };

  /**
   * \brief This is a specialization of the templated class for p = 0.
   *        This class represents the space \f$S^{0}_{0}\f$
   */
  template <>
  class ContinuousSpace<0> : public AbstractBEMSpace{
  public:
    static double b1(double t) {
      assert(IsWithinParameterRange(t));
      return 1.;
    }
    // Local to Global Map
    static unsigned int map(unsigned int q_, unsigned int n, unsigned int N) {
      assert(q_<=1 && n<=N);
      return 1;
    }
    // Constructor
    ContinuousSpace() {
      q = 1;
      referenceshapefunctions.push_back(b1);
      locglobmap = map;
    }
  };

  /**
   * \brief This is a specialization of the templated class for p = 1.
   *        This class represents the space \f$S^{0}_{1}\f$
   */
  template <>
  class ContinuousSpace<1> : public AbstractBEMSpace{
  public:
    // Basis function 1
    static double b1(double t) {
      assert(IsWithinParameterRange(t));
      return 0.5*(t+1);
    }
    // Basis function 2
    static double b2(double t) {
      assert(IsWithinParameterRange(t));
      return 0.5*(1-t);
    }
    // Local to Global Map
    static unsigned int map(unsigned int q_, unsigned int n, unsigned int N) {
      assert(q_<=2 && n<=N);
      if (q_==1)
        return n;
      else
        return ((n-1)%N==0) ? N : (n-1)%N==0;
    }
    // Constructor
    ContinuousSpace() {
      q = 2;
      referenceshapefunctions.push_back(b1);
      referenceshapefunctions.push_back(b2);
      locglobmap = map;
    }
  };

  /**
   * \brief This is a specialization of the templated class for p = 2.
   *        This class represents the space \f$S^{0}_{2}\f$
   */
  template <>
  class ContinuousSpace<2> : public AbstractBEMSpace{
  public:
    // Basis function 1
    static double b1(double t) {
      assert(IsWithinParameterRange(t));
      return 0.5*(t+1);
    }
    // Basis function 2
    static double b2(double t) {
      assert(IsWithinParameterRange(t));
      return 0.5*(1-t);
    }
    // Basis function 3
    static double b3(double t) {
      assert(IsWithinParameterRange(t));
      return 1-t*t;
    }
    static unsigned int map(unsigned int q_, unsigned int n, unsigned int N) {
      assert(q_<=3 && n<=N);
      if (q_==1)
        return n;
      else if (q_==2)
        return ((n-1)%N==0) ? N : (n-1)%N==0;
      else
        return N+n;
    }
    // Constructor
    ContinuousSpace() {
      q = 3;
      referenceshapefunctions.push_back(b1);
      referenceshapefunctions.push_back(b2);
      referenceshapefunctions.push_back(b3);
      locglobmap = map;
    }
  };
} // namespace parametricbem2d

#endif//CONTSPACEHPP
