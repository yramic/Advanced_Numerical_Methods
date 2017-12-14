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

    const double* gauss_pts_p = getGaussPoints(p);
    const double* gauss_wht_p = getGaussWeights(p);

    for(int i=1; i<=p; ++i) {

        for(int j=1; j<=p; ++j) {

            A(i-1,j-1) = 2.*std::sqrt(M_PI)*std::tgamma(1.+j) / ((3.+2.*i+2.*j)*std::tgamma(3./2.+j));
        }

        for(int k=0; k<p; ++k) {

            double tk = 0.5 * (gauss_pts_p[k] + 1.);
            double wk = 0.5 *  gauss_wht_p[k];

            b(i-1) += wk * std::pow(tk,i) * y(tk);
        }
    }

    Eigen::VectorXd x = A.colPivHouseholderQr().solve(b);

    size_t N = std::round(1./tau);
    Eigen::VectorXd grid = Eigen::VectorXd::LinSpaced(N+1, 0., 1.);
    Eigen::VectorXd u    = Eigen::VectorXd::Zero(N+1);

    for(int i=0; i<N+1; ++i) {
        for(int j=1; j<=p; ++j) {
            u(i) += x(j-1) * std::pow(grid(i),j);
        }
    }

    return grid;
}


Eigen::MatrixXcd toeplitz_triangular(const Eigen::VectorXcd& c)
{
    size_t n = c.size();
    Eigen::MatrixXcd T = Eigen::MatrixXcd::Zero(n, n);
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
Eigen::VectorXd cq_ieul_abel(const FUNC& y, size_t N)
{
    Eigen::VectorXcd w = Eigen::VectorXd::Zero(N+1);
    double tau = 1./N;

    int p = 8; // order of quadrature rule
    const double* gauss_pts_p = getGaussPoints(p);
    const double* gauss_wht_p = getGaussWeights(p);

    for(int i=0; i<p; ++i) {

        // integrate on semi-circumference centered in (1,0) with unitary radius:
        std::complex<double> ti = 1. + std::exp( std::complex<double>(0.,-0.5*M_PI*gauss_pts_p[i]) ); // M\_PI/2 - 0.5*(gauss\_pts\_p[i]+1.)*M\_PI
                     double  wi = 0.5 * M_PI * gauss_wht_p[i]; // change of integration domain to semi-circumference

        for(int j=0; j<N+1; ++j) {
            w(j) += wi / (std::sqrt(ti)*std::pow(1.-tau*ti,j+1));
        }

        // integrate on segment from (+1,-1) to (+1,+1):
        ti = std::complex<double>(1.,gauss_pts_p[i]);
        wi = gauss_wht_p[i];

        for(int j=0; j<N+1; ++j) {
            w(j) += wi / (std::sqrt(ti)*std::pow(1.-tau*ti,j+1));
        }
    }

    w *= std::tgamma(0.5)*tau / std::complex<double>(0.,2.*M_PI);

    // Solve the convolution quadrature:

    Eigen::VectorXd  grid = Eigen::VectorXd::LinSpaced(N+1,0.,1.);
    Eigen::VectorXcd y_N(N+1);
    for(int i=0; i<N+1; ++i) {
        y_N(i) = std::complex<double>(y(grid(i)),0.);
    }

    Eigen::MatrixXcd T = toeplitz_triangular(w);
    Eigen::VectorXcd u = T.triangularView<Eigen::Lower>().solve(y_N);

    return u.real();
}


/* @brief Find the unknown function u in the Abel integral equation
 * using convolution quadrature (BDF-2)
 * \param y Template function for the right-hand side
 * \param N Number of discretization steps
 * \\return Values of u from convolution quadrature
 */
template<typename FUNC>
Eigen::VectorXd cq_bdf2_abel(const FUNC& y, size_t N)
{
    Eigen::VectorXcd w = Eigen::VectorXd::Zero(N+1);
    double tau = 1./N;

    int p = 8; // order of quadrature rule
    const double* gauss_pts_p = getGaussPoints(p);
    const double* gauss_wht_p = getGaussWeights(p);

    for(int i=0; i<p; ++i) {

        // integrate on semi-circumference centered in (1,0) with unitary radius:
        std::complex<double> ti = 1. + std::exp( std::complex<double>(0.,-0.5*M_PI*gauss_pts_p[i]) ); // M\_PI/2 - 0.5*(gauss\_pts\_p[i]+1.)*M\_PI
                     double  wi = 0.5 * M_PI * gauss_wht_p[i]; // change of integration domain to semi-circumference

        for(int j=0; j<N+1; ++j) {
            w(j) += wi / (std::sqrt(ti)*std::sqrt(1.+tau*ti)) * (1./std::pow(2.-std::sqrt(1.+2.*tau*ti),j+1) - 1./std::pow(2.+std::sqrt(1.+2.*tau*ti),j+1));
        }

        // integrate on segment from (+1,-1) to (+1,+1):
        ti = std::complex<double>(1.,gauss_pts_p[i]);
        wi = gauss_wht_p[i];

        for(int j=0; j<N+1; ++j) {
            w(j) += wi / (std::sqrt(ti)*std::sqrt(1.+tau*ti)) * (1./std::pow(2.-std::sqrt(1.+2.*tau*ti),j+1) - 1./std::pow(2.+std::sqrt(1.+2.*tau*ti),j+1));
        }
    }

    w *= std::tgamma(0.5)*tau / std::complex<double>(0.,2.*M_PI);

    // Solve the convolution quadrature:

    Eigen::VectorXd  grid = Eigen::VectorXd::LinSpaced(N+1,0.,1.);
    Eigen::VectorXcd y_N(N+1);
    for(int i=0; i<N+1; ++i) {
        y_N(i) = std::complex<double>(y(grid(i)),0.);
    }

    Eigen::MatrixXcd T = toeplitz_triangular(w);
    Eigen::VectorXcd u = T.triangularView<Eigen::Lower>().solve(y_N);

    return u.real();
}


int main() {
    {
        auto y = [](double t) { return t; };

        double tau = 0.01;
        size_t N = std::round(1./tau);
        Eigen::VectorXd grid = Eigen::VectorXd::LinSpaced(N+1,0.,1.);
        Eigen::VectorXd u_ex(N+1);
        for(int i=0; i<N+1; ++i) {
            u_ex(i) = 8./M_PI*std::sqrt(grid(i));
        }

        std::cout << "Problem 3.1.g" << std::endl;
        for(int p=2; p<=32; p*=2) {
            Eigen::VectorXd u_app = poly_spec_abel(y, p, tau);
            Eigen::VectorXd diff  = u_ex - u_app;
            double err_max  = diff.cwiseAbs().maxCoeff();
            std::cout <<   "p = " << p << std::setw(15)
                      << "Max = "
                      << std::scientific << std::setprecision(3)
                      << err_max << std::endl;
        }
    }

    {
        auto y = [](double t) { return t; };

        std::cout << "Problem 3.1.l"  << std::endl;
        std::cout << "Implicit Euler" << std::endl;
        for(int N=10; N<=1280; N*=2) {

            Eigen::VectorXd grid = Eigen::VectorXd::LinSpaced(N+1,0.,1.);
            Eigen::VectorXd u_ex(N+1);
            for(int i=0; i<N+1; ++i) {
                u_ex(i) = 8./M_PI*std::sqrt(grid(i));
            }

            Eigen::VectorXd u_app = cq_ieul_abel(y, N);
            Eigen::VectorXd diff  = u_ex - u_app;
            double err_max  = diff.cwiseAbs().maxCoeff();
            std::cout <<   "N = " << N << std::setw(15)
                      << "Max = "
                      << std::scientific << std::setprecision(3)
                      << err_max << std::endl;
        }

        std::cout << "BDF-2" << std::endl;
        for(int N=10; N<=1280; N*=2) {

            Eigen::VectorXd grid = Eigen::VectorXd::LinSpaced(N+1,0.,1.);
            Eigen::VectorXd u_ex(N+1);
            for(int i=0; i<N+1; ++i) {
                u_ex(i) = 8./M_PI*std::sqrt(grid(i));
            }

            Eigen::VectorXd u_app = cq_bdf2_abel(y, N);
            Eigen::VectorXd diff  = u_ex - u_app;
            double err_max  = diff.cwiseAbs().maxCoeff();
            std::cout <<   "N = " << N << std::setw(15)
                      << "Max = "
                      << std::scientific << std::setprecision(3)
                      << err_max << std::endl;
        }
    }
}
