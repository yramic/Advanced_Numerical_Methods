/**
 * \file discontinuous_space.hpp
 * \brief This file declares a templated class inherited from AbstractBEMSpace
 *        and represents discontinuous spaces of the form \f$S^{-1}_{p}\f$
 *        as defined in \f$\eqref{eq:Qp}\f$, using full template specialization.
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
 *        defined in \f$\eqref{eq:Qp}\f$. The class is implemented through full
 *        template specialization for different values of p.
 */
template <unsigned int p> class DiscontinuousSpace : public AbstractBEMSpace {
public:
  DiscontinuousSpace() {
    std::cout << "Error! No specialization defined for p = " << p << std::endl;
  }
};

/**
 * \brief This is a specialization of the templated class for p = 0.
 *        This class represents the space \f$S^{-1}_{0}\f$
 */
template <> class DiscontinuousSpace<0> : public AbstractBEMSpace {
public:
  // Local to Global Map
  unsigned int LocGlobMap(unsigned int q, unsigned int n,
                          unsigned int N) const {
    // Asserting the index of local shape function and the panel number are
    // within limits
    assert(q <= q_ && n <= N);
    return n;
  }
  // Space Dimensions as defined in \f$\ref{T:thm:dimbe}\f$
  unsigned int getSpaceDim(unsigned int numpanels) const {
    return numpanels * q_;
  }
  // Constructor
  DiscontinuousSpace() {
    // Number of reference shape functions for the space
    q_ = 1;
    // Reference shape function 1, defined using a lambda expression
    BasisFunctionType b1 = [&](double t) { return 1.; };
    // Adding the reference shape functions to the vector
    referenceshapefunctions_.push_back(b1);
    // Reference shape function 1 derivative, defined using a lambda expression
    BasisFunctionType b1dot = [&](double t) { return 0.; };
    // Adding the reference shape function derivatives to the vector
    referenceshapefunctiondots_.push_back(b1dot);
  }
};

/**
 * \brief This is a specialization of the templated class for p = 1.
 *        This class represents the space \f$S^{-1}_{1}\f$
 */
template <> class DiscontinuousSpace<1> : public AbstractBEMSpace {
public:
  // Local to Global Map
  unsigned int LocGlobMap(unsigned int q, unsigned int n,
                          unsigned int N) const {
    // Asserting the index of local shape function and the panel number are
    // within limits
    assert(q <= q_ && n <= N);
    if (q == 1)
      return n;
    else
      return N + n;
  }
  // Space Dimensions as defined in \f$\ref{T:thm:dimbe}\f$
  unsigned int getSpaceDim(unsigned int numpanels) const {
    return numpanels * q_;
  }
  // Constructor
  DiscontinuousSpace() {
    // Number of reference shape functions for the space
    q_ = 2;
    // Reference shape function 1, defined using a lambda expression
    BasisFunctionType b1 = [&](double t) { return 0.5; };
    // Reference shape function 2, defined using a lambda expression
    BasisFunctionType b2 = [&](double t) { return 0.5*t; };
    // Adding the reference shape functions to the vector
    referenceshapefunctions_.push_back(b1);
    referenceshapefunctions_.push_back(b2);
    // Reference shape function 1 derivative, defined using a lambda expression
    BasisFunctionType b1dot = [&](double t) { return 0.; };
    // Reference shape function 2 derivative, defined using a lambda expression
    BasisFunctionType b2dot = [&](double t) { return 0.5; };
    // Adding the reference shape function derivatives to the vector
    referenceshapefunctiondots_.push_back(b1dot);
    referenceshapefunctiondots_.push_back(b2dot);
  }
};
} // namespace parametricbem2d

#endif // DISCONTSPACEHPP
