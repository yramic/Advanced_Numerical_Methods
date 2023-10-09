/***********************************************************************
 *                                                                     *
 * Code for Course "Advanced Numerical Methods for CSE"                *
 * (Prof. Dr. R. Hiptmair)                                             *
 * Author: Daniele Casati                                              *
 * Date: 11/2017                                                       *
 * (C) Seminar for Applied Mathematics, ETH Zurich                     *
 * This code can be freely used for non-commercial purposes as long    *
 * as this header is left intact.                                      *
 ***********************************************************************/
#ifndef KERNEL_HPP
#define KERNEL_HPP

/**
 * \brief Kernel functor \f$\frac{C}{|x-y|}\f$ if \f$x != y\f$, else 0
 */
class Kernel {
 public:
  /**
   * \brief Default constructor
   */
  Kernel() : num_(0.) {}

  /**
   * \brief Constructor
   * \param C Coefficient of the kernel function
   */
  Kernel(double num) : num_(num) {}

  /**
   * \brief Functor
   * \param x x-coordinate of grid point
   * \param y y-coordinate of grid point
   */
  virtual double operator()(double x, double y) = 0;

 protected:
  double num_;  //!< coefficient
};

/*!
 * \brief Kernel functor \f$C\log{\left|x-y\right|}\f$
 */
class KernelLog : public Kernel {
  using Kernel::Kernel;  // C++11 inheritance of constructors

 public:
  /**
   * \brief Functor
   * \param x x-coordinate of grid point
   * \param y y-coordinate of grid point
   */
  double operator()(double x, double y);
};

/*!
 * \brief Kernel functor \f$x \cdot y\f$:
 * low-rank approximation should be exact!
 */
class KernelPolynomial : public Kernel {
  using Kernel::Kernel;  // C++11 inheritance of constructors

 public:
  /**
   * \brief Functor
   * \param x x-coordinate of grid point
   * \param y y-coordinate of grid point
   */
  double operator()(double x, double y);
};

/*!
 * \brief Kernel functor \f$\frac{C}{\left|x-y\right|}\f$
 */
class KernelInvDistance : public Kernel {
  using Kernel::Kernel;  // C++11 inheritance of constructors

 public:
  /**
   * \brief Functor
   * \param x x-coordinate of grid point
   * \param y y-coordinate of grid point
   */
  double operator()(double x, double y);
};

/*!
 * \brief Kernel functor \f$\frac{\cos(C\left|x-y\right|)}{\left|x-y\right|}\f$
 */
class KernelCosine : public Kernel {
  using Kernel::Kernel;  // C++11 inheritance of constructors

 public:
  /**
   * \brief Functor
   * \param x x-coordinate of grid point
   * \param y y-coordinate of grid point
   */
  double operator()(double x, double y);
};

#endif  // KERNEL_HPP
