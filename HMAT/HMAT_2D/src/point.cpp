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
#include "../include/point.hpp"

// set X coordinate of the point
void Point::setX( double x) {
    x_ = x;
}

// set Y coordinate of the point
void Point::setY( double y) {
    y_ = y;
}

// set ID of the point
void Point::setId( double id) {
    id_ = id;
}

// set value of the point (not used)
void Point::setV( double v) {
    v_ = v;
}
