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
    //virtual ~Point();

    /**
    * \brief Getters
    */
    // return x coordinate of this point
    double getX() const {
        return x_;
    }
    // return y coordinate of this point
    double getY() const {
        return y_;
    }
    // return id of this point
    double getId() const {
        return id_;
    }
    // return value of this point
    double getV() const {
        return v_;
    }


    /**
    * \brief Setters
    */
    // set x coordinate
    void setX( double x) {
        x_ = x;
    };
    // set y coordinate
    void setY( double y) {
        y_ = y;
    };
    // set id
    void setId( double id) {
        id_ = id;
    };
    // set value
    void setV( double v) {
        v_ = v;
    };

private:

    double x_,y_;
    unsigned id_;
    double v_;
};

#endif // POINT_H
