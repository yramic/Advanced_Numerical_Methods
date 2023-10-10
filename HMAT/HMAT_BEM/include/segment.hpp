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
#ifndef SEGMENT_HPP
#define SEGMENT_HPP

#include <Eigen/Dense>
#include <vector>

/*!
* \brief Class for point attributes
*/
class Segment {
 public:
  /*!
    * \brief Default Constructor
    */
  Segment()
      : a_(Eigen::Vector2d::Zero()), b_(Eigen::Vector2d::Zero()), id_(0) {}
  /*!
    * \brief Return mean coordinates of this segment
    */
  Eigen::Vector2d getMean() const { return (a_ + b_) / 2.; }
  /*!
    * \brief Return first  extreme coordinates of this segment
    */
  Eigen::Vector2d getA() const { return a_; }
  /*!
    * \brief Return second extreme coordinates of this segment
    */
  Eigen::Vector2d getB() const { return b_; }
  /*!
    * \brief Return id of this point
    */
  unsigned getId() const { return id_; }
  /*!
    * \brief Set x coordinate
    */
  void setA(const Eigen::Vector2d& a);
  /*!
    * \brief Set y coordinate
    */
  void setB(const Eigen::Vector2d& b);
  /*!
    * \brief Set id
    */
  void setId(unsigned id);

 private:
  Eigen::Vector2d a_, b_;  //!< a, b extremes of the segment
  unsigned id_;            //!< id of the segment
};

#endif  // SEGMENT_HPP
