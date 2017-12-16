//// 
//// Copyright (C) 2017 SAM (D-MATH) @ ETH Zurich
//// Author(s): curzuato < > 
//// Contributors:  dcasati 
//// This file is part of the AdvNumCSE repository.
////
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <Eigen/SparseLU>
#include <cmath>
#include <iomanip>
#include <iostream>

using namespace Eigen;
using namespace std;


/* @brief Compute the convolution weights for Laplace transform F
 * \param F Template function for the Laplace transform
 * \param n Number of convolution weights, minus 1
 * \param p Order of quadrature rule
 * \param r Radius of circumference of integration (default = 1.0E-7)
 * \\return Values of u at final time t = 1
 */
template<typename FFUNC>
VectorXd conv_wght_bdf2(const FFUNC& F, size_t n, int p, double r=1.0E-7)
{
    // TODO: Compute the convolution weights for Laplace transform F
}


/* @brief Find the unknown function u at final time t = 1 in the evolution problem
 * using Galerkin discretization and convolution quadrature (BDF-2)
 * \param g Template function for the right-hand side
 * \param M Number of discretization intervals in time
 * \param N Number of discretization intervals in space
 * \param p Order of quadrature rule
 * \\return Values of u at final time t = 1
 */
template<typename FUNC>
VectorXd solveABC(const FUNC& g, size_t M, size_t N, int p)
{
    // TODO: Find the unknown function u at final time t = 1 in the evolution problem
}


int main() {
    // TODO: Tabulate the H1-error of the Galerkin discretization + convolution quadrature
}
