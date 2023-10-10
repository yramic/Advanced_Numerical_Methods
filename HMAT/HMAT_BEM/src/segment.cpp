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
#include "../include/segment.hpp"

// set X coordinate of the point
void Segment::setA(const Eigen::Vector2d& a)
{
    a_ = a;
}

// set Y coordinate of the point
void Segment::setB(const Eigen::Vector2d& b)
{
    b_ = b;
}

// set ID of the point
void Segment::setId(unsigned id)
{
    id_ = id;
}
