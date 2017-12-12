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


/* @brief Find the unknown function u in the Abel integral equation
 * using Galerkin discretization with a polynomial basis.
 * \param y Template function for the right-hand side
 * \param p Maximum degree of the polynomial basis and
 * order of the quadrature rule to compute the righ-hand side
 * \param tau Meshwidth of the grid where to compute the values of u
 * \\return Values of u on a grid in [0,1] with meshwidth tau
 */
template<typename FUNC>
Eigen::VectorXd poly_spec_abel(const FUNC& y, size_t p, double tau)
{
    // TODO: Find the unknown function u in the Abel integral equation
}


int main() {
    // TODO: Tabulate the max error of the Galerkin approximation scheme
}
