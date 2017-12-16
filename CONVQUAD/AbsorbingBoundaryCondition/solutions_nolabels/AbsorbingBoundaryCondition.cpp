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


/* @brief Find the unknown function u at final time t = 1 in the evolution problem
 * using Galerkin discretization and convolution quadrature (BDF-2)
 * \param phi Template function for the right-hand side
 * \param M Number of discretization intervals in time
 * \param N Number of discretization intervals in space
 * \param p Order of quadrature rule
 * \\return Values of u at final time t = 1
 */
template<typename FUNC>
VectorXd solveABC(const FUNC& phi, size_t M, size_t N, int p)
{
    double h = 1./N;
    VectorXd grid = VectorXd::LinSpaced(N+1,0.,1.);

    SparseMatrix<double> A(N+1,N+1); A.reserve(3*N+1); // 3*(N+1) - 2
    A.insert(0,0) = 1./h - 1./(h*h*pow(M_PI,3))*(2. + (M_PI*M_PI*h*h-2.)*cos(M_PI*h) - 2.*M_PI*h*sin(M_PI*h));
    A.insert(1,0) =-1./h + 1./(h*h*pow(M_PI,3))*(2. - 2.*cos(M_PI*h) - M_PI*h*sin(M_PI*h));
    for(int i=1; i<N; ++i) {
        A.insert(i-1,i) =-1./h + 1./(h*h*pow(M_PI,3))*(2.*(cos(M_PI*grid(i-1)) - cos(M_PI*grid(i))) - M_PI*h*(sin(M_PI*grid(i-1)) + sin(M_PI*grid(i))));
        A.insert(i,  i) = 2./h + 1./(h*h*pow(M_PI,3))*((M_PI*M_PI*h*h-4.)*(cos(M_PI*grid(i-1)) - cos(M_PI*grid(i+1))) + 4.*M_PI*h*(sin(M_PI*grid(i-1)) + sin(M_PI*grid(i+1))));
        A.insert(i+1,i) =-1./h + 1./(h*h*pow(M_PI,3))*(2.*(cos(M_PI*grid(i)) - cos(M_PI*grid(i+1))) - M_PI*h*(sin(M_PI*grid(i)) + sin(M_PI*grid(i+1))));
    }
    A.insert(N-1,N) =-1./h + 1./(h*h*pow(M_PI,3))*(2.*cos(M_PI*(1.-h)) - 2. - M_PI*h*sin(M_PI*(1.-h)));
    A.insert(N,  N) = 1./h + 1./(h*h*pow(M_PI,3))*((M_PI*M_PI*h*h-2.)*cos(M_PI*(1.-h)) - 2. + 2.*M_PI*h*sin(M_PI*(1.-h)));

    double tau = 1./M;
    VectorXcd w = VectorXd::Zero(M+1);

    for(int i=0; i<p; ++i) {

        // integrate on circumference centered in (0,0) with radius $~ 10^{-7}$:
        complex<double> ti = 1e-7*complex<double>(cos(2.*M_PI*i/p),sin(2.*M_PI*i/p));

        for(int l=0; l<M+1; ++l) {
            w(l) += 1./sqrt(1.+tau*ti) * (1./pow(2.-sqrt(1.+2.*tau*ti),l+1) - 1./pow(2.+sqrt(1.+2.*tau*ti),l+1)) * log(ti)/(ti*ti+1.);
        }
    }

    w *= tau*1e-7/complex<double>(0.,p); // coefficient of weight formula, tau/complex<double>(0.,2.*M\_PI), times quadrature weight, 2.*M\_PI*1e-7/p

    SparseMatrix<complex<double> > Aw = A.cast<complex<double> >();
    Aw.coeffRef(N,N) += w(0);
    SparseLU<SparseMatrix<complex<double> > > solver;
    solver.compute(Aw);

    MatrixXcd u = MatrixXcd::Zero(N+1,M+1); // first column is at time t = 0
    for(int i=1; i<M+1; ++i) {

        complex<double> rhs_N = 0.;
        for(int l=0; l<i; ++l) {
            rhs_N += w(i-l) * u(N,l);
        }

        VectorXcd phi_i(N+1);
        for(int j=0; j<N+1; ++j) {
            phi_i(j) = (complex<double>)phi(grid(j),i*tau);
        }

        VectorXcd rhs = phi_i; rhs(N) -= rhs_N;

        VectorXcd tmp = solver.solve(rhs);
        u.col(i) = tmp;
    }

    return u.col(M).real();
}


int main() {
    auto phi = [](double x, double t) { return sin(M_PI*x)*sin(M_PI*t); };
    VectorXd u_ref = solveABC(phi, 2560, 2560, 20);

    cout << "Problem 3.2.d" << endl;
    for(int N=10; N<=1280; N*=2) {
        VectorXd u_tmp = solveABC(phi, 2560, N, 20);
        double error = 0.;
        double ratio = 2560/N;
        for(int i=1; i<=2560; ++i) {
            int j = ceil(i/ratio);
            error += pow((u_ref(i) - u_ref(i-1))/2560. - (u_tmp(j) - u_tmp(j-1))/(double)(N-1), 2) / 2560.; // denominator is h of refined grid
        }
        cout << "N = " << N << setw(15)
             << "H1-Error = "
             << scientific << setprecision(3)
             << error << endl;
    }

    for(int M=10; M<=1280; M*=2) {
        VectorXd u_tmp = solveABC(phi, M, 2560, 20);
        double error = 0.;
        for(int i=1; i<=2560; ++i) {
            error += pow((u_ref(i) - u_ref(i-1))/2560. - (u_tmp(i) - u_tmp(i-1))/2560., 2) / 2560.; // denominator is h of refined grid
        }
        cout << "M = " << M << setw(15)
             << "H1-Error = "
             << scientific << setprecision(10)
             << error << endl;
    }
}
