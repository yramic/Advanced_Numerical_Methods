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
#ifndef CHEBY_HPP
#define CHEBY_HPP

#include <Eigen/Dense>

/**
 * \brief Compute Chebyshew nodes "tk_" and weights "wk_"
 * on domain [xl,xr] for a certain polynomial degree
 */
class Cheby {
 public:
  /*!
   * \brief Constructor
   * \param xl left coordinate of the interpolating domain
   * \param xr right coordinate of the interpolating domain
   */
  Cheby(double xl, double xr, unsigned deg);
  /*!
   * \brief Return Chebyshew nodes on domain [xl,xr]
   */
  Eigen::VectorXd getNodes() const { return tk_; }
  /*!
   * \brief Return weights of Lagrange polynomial
   */
  Eigen::VectorXd getWghts() const { return wk_; }
  /*!
   * \brief Compute Chebyshew nodes on domain [xl,xr]
   */
  void setNodes();
  /*!
   * \brief Compute weights of Lagrange polynomial
   */
  void setWghts();

 private:
  double xl_;  //!< left  boundary of domain on which we compute the Chebyshew
               //!< nodes
  double xr_;  //!< right boundary of domain on which we compute the Chebyshew
               //!< nodes
  unsigned deg_;        //!< degree of Lagrange polynomial
  Eigen::VectorXd tk_;  //!< Chebyshew nodes
  Eigen::VectorXd wk_;  //!< weights of Lagrange polynomial
};

#endif  // CHEBY_HPP
