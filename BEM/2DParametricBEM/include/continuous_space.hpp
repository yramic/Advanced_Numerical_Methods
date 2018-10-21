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
   * \brief This is a specialization of the templated class for p = 1.
   *        This class represents the space \f$S^{0}_{1}\f$
   */
  template <>
  class ContinuousSpace<1> : public AbstractBEMSpace{
  public:
    // Local to Global Map
    unsigned int LocGlobMap(unsigned int q, unsigned int n, unsigned int N) const{
      assert(q<=q_ && n<=N);
      if (q==2)
        return n;
      else
        return (n%N==0) ? 1 : (n+1);
    }
    // Space Dimensions
    unsigned int getSpaceDim(unsigned int numpanels) const{
      return numpanels*(q_-1);
    }
    // Constructor
    ContinuousSpace() {
      q_ = 2;
      // Basis function 1
      BasisFunctionType b1 = [&](double t) {
        return 0.5*(t+1);
      };
      // Basis function 2
      BasisFunctionType b2 = [&](double t) {
        return 0.5*(1-t);
      };
      referenceshapefunctions_.push_back(b1);
      referenceshapefunctions_.push_back(b2);
    }
  };

  /**
   * \brief This is a specialization of the templated class for p = 2.
   *        This class represents the space \f$S^{0}_{2}\f$
   */
  template <>
  class ContinuousSpace<2> : public AbstractBEMSpace{
  public:
    unsigned int LocGlobMap(unsigned int q, unsigned int n, unsigned int N) const{
      assert(q<=q_ && n<=N);
      if (q==2)
        return n;
      else if (q==1)
        return (n%N==0) ? 1 : (n+1);
      else
        return N+n;
    }
    // Space Dimensions
    unsigned int getSpaceDim(unsigned int numpanels) const{
      return numpanels*(q_-1);
    }
    // Constructor
    ContinuousSpace() {
      q_ = 3;
      // Basis function 1
      BasisFunctionType b1 = [&](double t) {
        return 0.5*(t+1);
      };
      // Basis function 2
      BasisFunctionType b2 = [&](double t) {
        return 0.5*(1-t);
      };
      // Basis function 3
      BasisFunctionType b3 = [&](double t) {
        return 1-t*t;
      };
      referenceshapefunctions_.push_back(b1);
      referenceshapefunctions_.push_back(b2);
      referenceshapefunctions_.push_back(b3);
    }
  };
} // namespace parametricbem2d

#endif//CONTSPACEHPP
