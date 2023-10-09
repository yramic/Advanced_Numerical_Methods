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
#ifndef POINT_HPP
#define POINT_HPP

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
  Point() : x_(0), y_(0), id_(0), v_(0) {}
  /*!
   * \brief Return x coordinate of this point
   */
  double getX() const { return x_; }
  /*!
   * \brief Return y coordinate of this point
   */
  double getY() const { return y_; }
  /*!
   * \brief Return id of this point
   */
  unsigned getId() const { return id_; }
  /*!
   * \brief Return value of this point (for debugging)
   */
  double getV() const { return v_; }
  /*!
   * \brief Set x coordinate
   */
  void setX(double x);
  /*!
   * \brief Set y coordinate
   */
  void setY(double y);
  /*!
   * \brief Set id
   */
  void setId(unsigned id);
  /*!
   * \brief Set value (for debugging)
   */
  void setV(double v);

 private:
  double x_, y_;  //!< x, y coordinate of the point
  unsigned id_;   //!< id of the point
  double v_;
};

#endif  // POINT_HPP
