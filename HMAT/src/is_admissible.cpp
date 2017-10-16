#include "../include/is_admissible.hpp"
#include <cmath>


double AdmissibilityD::get_max(double xl, double xr, double yl, double yr)
{
    if(std::abs(xr-xl) > std::abs(yr-yl))
        return std::abs(xr-xl);
    else
        return std::abs(yr-yl);
}


double AdmissibilityD::get_min(double xl, double xr, double yl, double yr)
{
    double dist = std::abs(xl-yl);
    double dist2 = std::abs(xl-yr);
    double dist3 = std::abs(xr-yl);
    double dist4 = std::abs(xr-yr);

    if(dist2 < dist)
        dist = dist2;
    if(dist3 < dist)
        dist = dist3;
    if(dist4 < dist)
        dist = dist4;

    return dist;
}


bool AdmissibilityD::is_admissible(double xl, double xr, double yl, double yr, double eta)
{
    return get_max(xl,xr,yl,yr) <= eta * get_min(xl,xr,yl,yr);
}

// find the biggest edge of each Bounding Box and compare it with the biggest node of the other Bounding Box
double AdmissibilityH::get_max(Node* a, Node* b)
{
    double xa,ya,xb,yb,maxA,maxB;
    xa = std::abs(a->getX2_b()-a->getX1_b());
    ya = std::abs(a->getY2_b()-a->getY1_b());
    if(xa > ya){
        maxA = xa;
    }
    else {
        maxA = ya;
    }
    xb = std::abs(b->getX2_b()-b->getX1_b());
    yb = std::abs(b->getY2_b()-b->getY1_b());
    if(xa > ya){
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
    double x1 = a->getX1_b();
    double x1b = a->getX2_b();
    double x2 = b->getX1_b();
    double x2b = b->getX2_b();
    double y1 = a->getY1_b();
    double y1b = a->getY2_b();
    double y2 = b->getY1_b();
    double y2b = b->getY2_b();
    bool left, right, bottom, top;
    left = x2b < x1;
    right = x1b < x2;
    bottom = y2b < y1;
    top = y1b < y2;
    if (top && left){
        return dist(x1, y1b, x2b, y2);
    }
    else if (left && bottom){
        return dist(x1, y1, x2b, y2b);
    }
    else if (bottom && right){
        return dist(x1b, y1, x2, y2b);
    }
    else if (right && top){
        return dist(x1b, y1b, x2, y2);
    }
    else if (left) {
        return (x1 - x2b);
    }
    else if (right) {
        return (x2 - x1b);
    }
    else if (bottom) {
        return (y1 - y2b);
    }
    else if (top) {
        return (y2 - y1b);
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
