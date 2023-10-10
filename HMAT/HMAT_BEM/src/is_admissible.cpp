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
#include "../include/is_admissible.hpp"
#include <cmath>

// find the biggest edge of each Bounding Box and compare it with the biggest node of the other Bounding Box
double AdmissibilityH::get_max(Node* a, Node* b)
{
    double xa,ya,xb,yb,maxA,maxB;
    xa = std::abs(a->getX1()-a->getX2());
    ya = std::abs(a->getY1()-a->getY2());
    if(xa > ya){
        maxA = xa;
    }
    else {
        maxA = ya;
    }
    xb = std::abs(b->getX1()-b->getX2());
    yb = std::abs(b->getY1()-b->getY2());
    if(xb > yb){
        maxB = xb;
    }
    else {
        maxB = yb;
    }
    if(maxA > maxB){
        return maxA;
    }
    else {
        return maxB;
    }
}

// return the distance between 2 points
double dist(double x, double y, double a, double b){
    return std::sqrt(std::pow(x - a, 2) + std::pow(y - b, 2));
}

// return the distance between 2 Boxes
double AdmissibilityH::get_min(Node* a, Node* b)
{
    double x1a = a->getX1();
    double x2a = a->getX2();
    double x1b = b->getX1();
    double x2b = b->getX2();
    double y1a = a->getY1();
    double y2a = a->getY2();
    double y1b = b->getY1();
    double y2b = b->getY2();
    bool left, right, bottom, top;
    left = x2b < x1a;
    right = x2a < x1b;
    bottom = y2b < y1a;
    top = y2a < y1b;
    if (top && left){
        return dist(x1a, y2a, x2b, y1b);
    }
    else if (left && bottom){
        return dist(x1a, y1a, x2b, y2b);
    }
    else if (bottom && right){
        return dist(x2a, y1a, x1b, y2b);
    }
    else if (right && top){
        return dist(x2a, y2a, x1b, y1b);
    }
    else if (left) {
        return (x1a - x2b);
    }
    else if (right) {
        return (x1b - x2a);
    }
    else if (bottom) {
        return (y1a - y2b);
    }
    else if (top) {
        return (y1b - y2a);
    }
    else {
        return 0;
    }
}

// return if two Bounding Boxes are admissible
bool AdmissibilityH::is_admissible(Node* x, Node* y, double eta)
{
    return get_max(x, y) <= eta * get_min(x, y);
}
