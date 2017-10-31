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
    Point():
        x_(0),y_(0),id_(0),v_(0)
    {}
    /*!
    * \brief Return x coordinate of this point
    */
    double getX() const {
        return x_;
    }
    /*!
    * \brief Return y coordinate of this point
    */
    double getY() const {
        return y_;
    }
    /*!
    * \brief Return id of this point
    */
    double getId() const {
        return id_;
    }
    /*!
    * \brief Return value of this point
    */
    double getV() const {
        return v_;
    }
    /*!
    * \brief Set x coordinate
    */
    void setX( double x) {
        x_ = x;
    }
    /*!
    * \brief Set y coordinate
    */
    void setY( double y) {
        y_ = y;
    }
    /*!
    * \brief Set id
    */
    void setId( double id) {
        id_ = id;
    }
    /*!
    * \brief Set value
    */
    void setV( double v) {
        v_ = v;
    }

private:
    double x_,y_;   //!< x, y coordinate of the point
    unsigned id_;   //!< id of the point
    double v_;
};

#endif // POINT_H
