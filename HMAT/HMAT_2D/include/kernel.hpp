/***********************************************************************
 *                                                                     *
 * Code for Course "Advanced Numerical Methods for CSE"                *
 * (Prof. Dr. R. Hiptmair)                                             *
 * Author: Ioannis Magkanaris                                          *
 * Date: 11/2017                                                       *
 * (C) Seminar for Applied Mathematics, ETH Zurich                     *
 * This code can be freely used for non-commercial purposes as long    *
 * as this header is left intact.                                      *
 ***********************************************************************/
#ifndef KERNEL_HPP
#define KERNEL_HPP

/*!
* \brief Kernel base class
*/
class Kernel {
 public:
  /*!
    * \brief Default Constructor
    */
  Kernel() : num_(0.) {}
  /*!
    * \brief Default Constructor
    */
  Kernel(double num) : num_(num) {}

  /*!
    * \brief Virtual Functor
    * \param x1 x-coordinate of first point
    * \param y1 y-coordinate of first point
    * \param x2 x-coordinate of second point
    * \param y2 y-coordinate of second point
    */
  virtual double operator()(double x1, double y1, double x2, double y2) = 0;

 protected:
  double num_;  //!< coefficient
};

/*!
* \brief Kernel functor \f$-\frac{1}{2 \pi}\log{\left|\vec{x}-\vec{y}\right|}\f$
*/
class KernelGalerkin : public Kernel {
  using Kernel::Kernel;  // C++11 inheritance of constructors

 public:
  /*!
    * \brief Functor for 2D
    * \param x1 x-coordinate of first point
    * \param y1 y-coordinate of first point
    * \param x2 x-coordinate of second point
    * \param y2 y-coordinate of second point
    */
  double operator()(double x1, double y1, double x2, double y2);
};

/*!
* \brief Kernel functor \f${x1}\times{x2}\times{y1}\times{y2}\f$
*/
class KernelPolynomial : public Kernel {
  using Kernel::Kernel;  // C++11 inheritance of constructors

 public:
  /*!
    * \brief Functor for 2D
    * \param x1 x-coordinate of first point
    * \param y1 y-coordinate of first point
    * \param x2 x-coordinate of second point
    * \param y2 y-coordinate of second point
    */
  double operator()(double x1, double y1, double x2, double y2);
};

/*!
* \brief Kernel functor 1
*/
class KernelConstant : public Kernel {
  using Kernel::Kernel;  // C++11 inheritance of constructors

 public:
  /*!
    * \brief Functor for 2D
    * \param x1 x-coordinate of first point
    * \param y1 y-coordinate of first point
    * \param x2 x-coordinate of second point
    * \param y2 y-coordinate of second point
    */
  double operator()(double x1, double y1, double x2, double y2);
};

/**
* \brief Kernel functor \f$\cos(|\vec{x}-\vec{y}|)\f$
*/
class KernelGlobalSmooth : public Kernel {
  using Kernel::Kernel;  // C++11 inheritance of constructors

 public:
  /**
    * \brief Functor for 2D
    * \param x1 x-coordinate of first point
    * \param y1 y-coordinate of first point
    * \param x2 x-coordinate of second point
    * \param y2 y-coordinate of second point
    */
  double operator()(double x1, double y1, double x2, double y2);
};

/**
* \brief Kernel functor \f$e^{-|\vec{x}-\vec{y}|^2}\f$
*/
class KernelGauss : public Kernel {
  using Kernel::Kernel;  // C++11 inheritance of constructors

 public:
  /**
    * \brief Functor for 2D
    * \param x1 x-coordinate of first point
    * \param y1 y-coordinate of first point
    * \param x2 x-coordinate of second point
    * \param y2 y-coordinate of second point
    */
  double operator()(double x1, double y1, double x2, double y2);
};

#endif  // KERNEL_HPP
