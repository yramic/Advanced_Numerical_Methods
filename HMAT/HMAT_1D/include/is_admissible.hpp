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
#ifndef IS_ADMISSIBLE_HPP
#define IS_ADMISSIBLE_HPP

/**
 * \brief Primitive functions used to check whether a cluster pair is admissible
 * (eta-admissibility)
 */

/*!
 * \brief Get max edge of the cluster
 * \param xl left x coordinate of cluster
 * \param xr right x coordinate of cluster
 * \param yl left y coordinate of cluster (corresponding to the x axis)
 * \param yr right x coordinate of cluster (corresponding to the x axis)
 */
double get_max(double xl, double xr, double yl, double yr);
/*!
 * \brief Get distance from diagonal
 * \param xl left x coordinate of cluster
 * \param xr right x coordinate of cluster
 * \param yl left y coordinate of cluster (corresponding to the x axis)
 * \param yr right x coordinate of cluster (corresponding to the x axis)
 */
double get_min(double xl, double xr, double yl, double yr);
/*!
 * \brief Check the admissibility of a bounding box
 * \param xl left x coordinate of cluster
 * \param xr right x coordinate of cluster
 * \param yl left y coordinate of cluster (corresponding to the x axis)
 * \param yr right x coordinate of cluster (corresponding to the x axis)
 * \param eta eta-admissibility constant
 */
bool is_admissible(double xl, double xr, double yl, double yr, double eta);

#endif // IS_ADMISSIBLE_HPP
