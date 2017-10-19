#ifndef GLOBAL_INTERPOLATION_APP_HPP
#define GLOBAL_INTERPOLATION_APP_HPP

#include "kernel.hpp"
#include "point.hpp"
#include <Eigen/Dense>

/**
* \brief Master class for global interpolation approximation
*/
class GlobalInterpolationApp
{
public:

    /**
    * \brief Constructor
    */
    //GlobalInterpolationApp(GlobalSmoothKernel gskernel, std::vector<Point> pp, int n);
    //GlobalInterpolationApp(ConstantKernel ckernel, std::vector<Point> pp, int n);
    //GlobalInterpolationApp(PolynomialKernel ckernel, std::vector<Point> pp, int n);
    //GlobalInterpolationApp(GaussKernel ckernel, std::vector<Point> pp, int n);
    GlobalInterpolationApp(Kernel* ckernel, std::vector<Point> pp, int n);
    // approximate matrix-vector multiplication
    Eigen::VectorXd mvProd(Eigen::VectorXd &c, unsigned deg);

private:
    GlobalSmoothKernel GSK_;    // Global Smooth Kernel
    ConstantKernel CK_;         // Constant Kernel
    PolynomialKernel PK_;       // Polynomial Kernel
    GaussKernel GK_;            // Gauss Kernel
    Kernel* K_;
    std::vector<Point> PPoints_;// Polygon Points
};

#endif // GLOBAL_INTERPOLATION_APP_HPP
