/**
 * \file discontinuous_space.hpp
 * \brief This file declares a templated class inherited from AbstractBEMSpace
 *        and represents discontinuous spaces of the form \f$S^{-1}_{p}\f$
 * (eq.1.4.22) using full template specialization.
 *
 * This File is a part of the 2D-Parametric BEM package
 */

#ifndef DISCONTSPACEHPP
#define DISCONTSPACEHPP

#include "abstract_bem_space.hpp"

#include <iostream>
#include <utility>
#include <vector>

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
class DiscontinuousSpace<0> : public AbstractBEMSpace {
public:
  // Local to Global Map
  unsigned int LocGlobMap(unsigned int q, unsigned int n,
                          unsigned int N) const {
    assert(q <= q_ && n <= N);
    return n;
  }
  // Space Dimensions
  unsigned int getSpaceDim(unsigned int numpanels) const {
    return numpanels * q_;
  }
  // Constructor
  DiscontinuousSpace() {
    q_ = 1;
    BasisFunctionType b1 = [&](double t) {
      return 1.;
    };
    referenceshapefunctions_.push_back(b1);
  }
};

/**
 * \brief This is a specialization of the templated class for p = 1.
 *        This class represents the space \f$S^{-1}_{1}\f$
 */
template <>
class DiscontinuousSpace<1> : public AbstractBEMSpace {
public:
  // Local to Global Map
  unsigned int LocGlobMap(unsigned int q, unsigned int n,
                          unsigned int N) const {
    assert(q <= q_ && n <= N);
    if (q == 1)
      return n;
    else
      return N + n;
  }
  // Space Dimensions
  unsigned int getSpaceDim(unsigned int numpanels) const {
    return numpanels * q_;
  }
  // Constructor
  DiscontinuousSpace() {
    q_ = 2;
    // Basis function 1
    BasisFunctionType b1 = [&](double t) { return 0.5; };
    // Basis function 2
    BasisFunctionType b2 = [&](double t) {
      return 0.5 * t;
    };
    referenceshapefunctions_.push_back(b1);
    referenceshapefunctions_.push_back(b2);
  }
};
} // namespace parametricbem2d

#endif // DISCONTSPACEHPP
