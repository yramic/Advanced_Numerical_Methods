//// 
//// Copyright (C) 2017 SAM (D-MATH) @ ETH Zurich
//// Author(s): curzuato < > 
//// Contributors:  dcasati 
//// This file is part of the AdvNumCSE repository.
////
#include <Eigen/Dense>
#include <cmath>
#include <iomanip>
#include <iostream>
extern "C" {
#include "../../../BEM/CppHilbert/Library/source/gaussQuadrature.h"
}

using namespace Eigen;
using namespace std;


/* @brief Find the unknown function u in the Abel integral equation
 * using Galerkin discretization with a polynomial basis.
 * \param y Template function for the right-hand side
 * \param p Maximum degree of the polynomial basis and
 * order of the quadrature rule to compute the righ-hand side
 * \param tau Meshwidth of the grid where to compute the values of u
 * \\return Values of u on a grid in [0,1] with meshwidth tau
 */
template<typename FUNC>
VectorXd poly_spec_abel(const FUNC& y, size_t p, double tau)
{
    // TODO: Find the unknown function u in the Abel integral equation with Galerkin discretization
}


MatrixXcd toeplitz_triangular(const VectorXcd& c)
{
    size_t n = c.size();
    MatrixXcd T = MatrixXcd::Zero(n, n);
    for(int i=0; i<n; ++i) {
        T.col(i).tail(n-i) = c.head(n-i);
    }
    return T;
}


/* @brief Find the unknown function u in the Abel integral equation
 * using convolution quadrature (implicit Euler)
 * \param y Template function for the right-hand side
 * \param N Number of discretization steps
 * \\return Values of u from convolution quadrature
 */
template<typename FUNC>
VectorXd cq_ieul_abel(const FUNC& y, size_t N)
{
    // TODO: Find the unknown function u in the Abel integral equation with convolution quadrature (implicit Euler)
}


/* @brief Find the unknown function u in the Abel integral equation
 * using convolution quadrature (BDF-2)
 * \param y Template function for the right-hand side
 * \param N Number of discretization steps
 * \\return Values of u from convolution quadrature
 */
template<typename FUNC>
VectorXd cq_bdf2_abel(const FUNC& y, size_t N)
{
    // TODO: Find the unknown function u in the Abel integral equation with convolution quadrature (BDF-2)
}


int main() {
    // TODO: Tabulate the max error of the Galerkin approximation scheme

    // TODO: Tabulate the max error of the convolution quadratures
}
