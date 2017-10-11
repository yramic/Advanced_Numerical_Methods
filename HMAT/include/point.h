#ifndef POINT_H
#define POINT_H

#include <Eigen/Dense>
#include <vector>

class Point
{
public:

    /**
    * \brief Constructors
    */
    // default constructor
    Point():
        x_(0),y_(0),id_(0),v_(0)
    {}


    // destructor
    virtual ~Point();

    /**
    * \brief Getters
    */
    // return x coordinate of this point
    double getx() const {
        return x_;
    }
    // return y coordinate of this point
    double gety() const {
        return y_;
    }
    // return id of this point
    double getid() const {
        return id_;
    }
    // return Value of this point
    double getV() const {
        return v_;
    }


    /**
    * \brief Setters
    */
    // set x coordinate
    void setx( double x) {
        x_ = x;
    };
    // set y coordinate
    void sety( double y) {
        y_ = y;
    };
    // set id
    void setid( double id) {
        id_ = id;
    };
    // set Value
    void setv( double v) {
        v_ = v;
    };

private:

    double x_,y_;
    unsigned id_;
    double v_;
};

#endif // POINT_H
