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
    Eigen::MatrixXd A = Eigen::MatrixXd::Zero(p,p);
    Eigen::VectorXd b = Eigen::VectorXd::Zero(p);

    const double* gauss_pts_pp = getGaussPoints(p*p);
    const double* gauss_wht_pp = getGaussWeights(p*p);

    const double* gauss_pts_p = getGaussPoints(p);
    const double* gauss_wht_p = getGaussWeights(p);

    for(int i=0; i<p; ++i) {

        for(int j=0; j<p; ++j) {

            for(int k=0; k<p*p; ++k) {

                double tk = 0.5 * (gauss_pts_pp[k] + 1.);
                double wk = 0.5 *  gauss_wht_pp[k];

                for(int l=0; l<p*p; ++l) {

                    double xl = 0.5 * std::sqrt(tk) * (gauss_pts_pp[k] + 1.);
                    double wl = 0.5 * std::sqrt(tk) *  gauss_wht_pp[k];

                    A(i,j) += wk*wl * 2. * std::pow(tk,i) * std::pow(tk - xl*xl,j);
                }
            }
        }

        for(int k=0; k<p; ++k) {

            double tk = 0.5 * (gauss_pts_p[k] + 1.);
            double wk = 0.5 *  gauss_wht_p[k];

            b(i) += wk * std::pow(tk,i) * y(tk);
        }
    }

    Eigen::VectorXd x = A.colPivHouseholderQr().solve(b);

    size_t N = std::round(1./tau) + 1;
    Eigen::VectorXd grid = Eigen::LinSpaced(N, 0., 1.);
    Eigen::VectorXd u    = Eigen::VectorXd::Zero(N);

    for(int i=0; i<N; ++i) {
        for(int j=0; j<p; ++j) {
            u(i) += x(j) * std::pow(grid(i),j);
        }
    }

    return grid;
}


int main() {
    double tau = 0.01;
    size_t N = std::round(1./tau) + 1;
    Eigen::VectorXd grid = Eigen::LinSpaced(N, 0., 1.);
    Eigen::VectorXd u_ex(N);
    for(int i=0; i<N; ++i) {
        u_ex(i) = 8./M_PI*std::sqrt(grid(i));
    }
    auto y = [](double t)  { return t; };

    std::cout << "Problem 3.1.g" << std::endl;
    for(int p=2; p<=10; ++p) {
        Eigen::VectorXd u_app = poly_spec_abel(y, p, tau);
        Eigen::VectorXd diff  = u_ex - u_app;
        double err_max  = diff.cwiseAbs().maxCoeff();
        std::cout <<   "p = " << p << std::setw(15)
                  << "Max = "
                  << std::scientific << std::setprecision(3)
                  << err_max << std::endl;
    }
}
