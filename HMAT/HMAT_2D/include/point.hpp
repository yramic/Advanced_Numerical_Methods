#ifndef POINT_H
#define POINT_H

#include <Eigen/Dense>
#include <vector>
/*!
* \brief Class for point attributes
*/
class Point
{
public:

    /*!
    * \brief Default Constructor
    */
    // default constructor
    Point():
        x_(0),y_(0),id_(0),v_(0)
    {}


    // destructor
    //virtual ~Point();

    /*!
    * \brief return x coordinate of this point
    */
    // return x coordinate of this point
    double getX() const {
        return x_;
    }
    /*!
    * \brief return y coordinate of this point
    */
    // return y coordinate of this point
    double getY() const {
        return y_;
    }
    /*!
    * \brief return id of this point
    */
    // return id of this point
    double getId() const {
        return id_;
    }
    /*!
    * \brief return value of this point
    */
    // return value of this point
    double getV() const {
        return v_;
    }


    /*!
    * \brief Set x coordinate
    */
    // set x coordinate
    void setX( double x) {
        x_ = x;
    }
    /*!
    * \brief Set y coordinate
    */
    // set y coordinate
    void setY( double y) {
        y_ = y;
    }
    /*!
    * \brief Set id
    */
    // set id
    void setId( double id) {
        id_ = id;
    }
    /*!
    * \brief Set value
    */
    // set value
    void setV( double v) {
        v_ = v;
    }

private:

    double x_,y_;
    unsigned id_;
    double v_;
};

#endif // POINT_H
