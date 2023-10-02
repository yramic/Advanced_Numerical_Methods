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
#ifndef POINT_H
#define POINT_H

#include <Eigen/Dense>
#include <vector>
/*!
 * \brief Class for point attributes
 */
class Point {
 public:
  /*!
   * \brief Default Constructor
   */
  Point() : x_(0), id_(0) {}
  /*!
   * \brief Return x coordinate of this point
   */
  double getX() const { return x_; }
  /*!
   * \brief Return id of this point
   */
  double getId() const { return id_; }
  /*!
   * \brief Set x coordinate
   */
  void setX(double x);
  /*!
   * \brief Set id
   */
  void setId(double id);

 private:
  double x_;     //!< x coordinate of the point
  unsigned id_;  //!< id of the point
};

#endif  // POINT_H
