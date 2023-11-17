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
#ifndef IS_ADMISSIBLE_HPP
#define IS_ADMISSIBLE_HPP

#include "node.hpp"

/*!
* \brief Primitive functions used to check whether a cluster is admissible
* (eta-admissibility)
*/
class Admissibility {
 public:
  /*!
     * \brief Return maximum base
     */
  virtual double get_max(Node* x, Node* y) = 0;
  /*!
     * \brief Return minimum base
     */
  virtual double get_min(Node* x, Node* y) = 0;
  /*!
     * \brief Return if the cluster is admissible
     */
  virtual bool is_admissible(Node* x, Node* y, double eta) = 0;
};

/*!
* \brief Class for 2D admissibility problems(admissibility between 2 clusters´ rectangles)
*/
class AdmissibilityH : public Admissibility {
 public:
  /*!
     * \brief Returns the biggest edge of the two cluster rectangles
     * \param x Node of the first cluster rectangle
     * \param y Node of the second cluster rectangle
     */
  double get_max(Node* x, Node* y);
  /*!
     * \brief Returns the distance of the two cluster rectangles
     * \param x Node of the first cluster rectangle
     * \param y Node of the second cluster rectangle
     */
  double get_min(Node* x, Node* y);
  /*!
     * \brief Returns if the two cluster rectangles are admissible
     * \param x Node of the first cluster rectangle
     * \param y Node of the second cluster rectangle
     * \param eta eta admissibility variable
     */
  bool is_admissible(Node* x, Node* y, double eta);
};

#endif  // IS_ADMISSIBLE_HPP
