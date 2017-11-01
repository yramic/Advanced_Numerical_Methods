#ifndef GLOBAL_INTERPOLATION_APP_HPP
#define GLOBAL_INTERPOLATION_APP_HPP

#include "kernel.hpp"
#include "point.hpp"
#include <Eigen/Dense>

/*!
* \brief Master class for global interpolation approximation
*/
class GlobalInterpolationApp
{
public:
    /*!
     * \brief Constructor for 2D problem
     * \param kernel Kernel used for the matrix multiplication
     * \param pp Vector of points in space
     * \param n Number of points
     */
    GlobalInterpolationApp(Kernel* kernel, std::vector<Point> pp, int n);
    /*!
     * \brief Approximate matrix-vector multiplication
     * \param c Vector c
     * \param deg Degree of itnerpolation
     */
    Eigen::VectorXd mvProd(Eigen::VectorXd &c, unsigned deg);

private:
    Kernel* K_;                     //!< Kernel pointer
    std::vector<Point> PPoints_;    //!< Polygon Points
};

#endif // GLOBAL_INTERPOLATION_APP_HPP
