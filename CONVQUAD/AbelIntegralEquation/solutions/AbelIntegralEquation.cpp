#include <Eigen/Dense>
#include <unsupported/Eigen/FFT>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <ctime>
#include "gauleg.hpp"
#include "utilities.hpp"

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
/* SAM_LISTING_BEGIN_0 */
template<typename FUNC>
VectorXd poly_spec_abel(const FUNC& y, size_t p, double tau)
{
    MatrixXd A = MatrixXd::Zero(p+1,p+1);
    VectorXd b = VectorXd::Zero(p+1);
    
    // generate Gauss-Legendre points and weights
    int ng = p;
    Eigen::RowVectorXd gauss_pts_p, gauss_wht_p;
    std::tie(gauss_pts_p,gauss_wht_p) = gauleg(0., 1., ng);
    
    // set-up the Galerkin matrix and rhs vector
    for(int i=0; i<=p; i++)
    {
        for(int j=0; j<=p; j++)
            A(i,j) = 2.*sqrt(M_PI)*tgamma(1+j) / ((3.+2.*i+2.*j)*tgamma(3./2.+j)); // tgamma(1+j) == j! if j is integer
    
        for(int k=0; k<ng; k++)
        {
            double tk = gauss_pts_p(k);
            double wk = gauss_wht_p(k);
            b(i) += wk * pow(tk,i) * y(tk);
        }
    }

    // linear system solve using QR decomposition
    VectorXd x = A.colPivHouseholderQr().solve(b);
    
    size_t N = round(1./tau);
    VectorXd grid = VectorXd::LinSpaced(N+1, 0., 1.);
    VectorXd u    = VectorXd::Zero(N+1);
    
    // generate solution at grid points
    for(int i=0; i<=N; i++) {
        for(int j=0; j<=p; j++) {
            u(i) += x(j) * pow(grid(i),j);
        }
    }

    return u;
}
/* SAM_LISTING_END_0 */


/* @brief Find the unknown function u in the Abel integral equation
 * using convolution quadrature (implicit Euler)
 * \param y Template function for the right-hand side
 * \param N Number of discretization steps
 * \\return Values of u from convolution quadrature
 */
/* SAM_LISTING_BEGIN_2 */
template<typename FUNC>
VectorXd cq_ieul_abel(const FUNC& y, size_t N)
{
    VectorXd w(N+1); w(0) = 1.;
    for(int l=1; l<N+1; ++l) {
        w(l) = w(l-1) * (l - 0.5) / l; // denominator is factorial
    }
    w *= sqrt(M_PI/N);

    // Solve the convolution quadrature:

    VectorXd grid = VectorXd::LinSpaced(N+1,0.,1.);
    VectorXd y_N(N+1);
    for(int i=0; i<N+1; ++i) {
        y_N(i) = y(grid(i));
    }
    
    MatrixXd T = toeplitz_triangular(w);
    VectorXd u = T.triangularView<Lower>().solve(y_N);
    
    return u;
}
/* SAM_LISTING_END_2 */


/* @brief Find the unknown function u in the Abel integral equation
 * using convolution quadrature (BDF-2)
 * \param y Template function for the right-hand side
 * \param N Number of discretization steps
 * \\return Values of u from convolution quadrature
 */
/* SAM_LISTING_BEGIN_3 */
template<typename FUNC>
VectorXd cq_bdf2_abel(const FUNC& y, size_t N)
{
    VectorXd w1(N+1); w1(0) = 1.;
    for(int l=1; l<N+1; ++l) {
        w1(l) = w1(l-1) * (l - 0.5) / l; // denominator is factorial
    }

    VectorXd w2 = w1;
    for(int l=1; l<N+1; ++l) {
        w2(l) /= pow(3,l);
    }

    VectorXd w = myconv(w1, w2).head(N+1).real();
    w *= sqrt(M_PI/N) * sqrt(2./3.);

    // Solve the convolution quadrature:

    VectorXd grid = VectorXd::LinSpaced(N+1,0.,1.);
    VectorXd y_N(N+1);
    for(int i=0; i<N+1; ++i) {
        y_N(i) = y(grid(i));
    }

    MatrixXd T = toeplitz_triangular(w);
    VectorXd u = T.triangularView<Lower>().solve(y_N);

    return u;
}
/* SAM_LISTING_END_3 */


int main() {
    /* SAM_LISTING_BEGIN_1 */
    {
        auto u = [](double t) { return 2./M_PI*sqrt(t); };
        auto y = [](double t) { return t; };
        
        double tau = 0.01;
        size_t N = round(1./tau);
        VectorXd grid = VectorXd::LinSpaced(N+1,0.,1.);
        VectorXd u_ex(N+1);
        for(int i=0; i<N+1; ++i) {
            u_ex(i) = u(grid(i));
        }
    
        cout << "\nSpectral Galerkin\n" << endl;
        for(int p=2; p<=10; ++p) {
            VectorXd u_app = poly_spec_abel(y, p, tau);
            VectorXd diff  = u_ex - u_app;
            double err_max = diff.cwiseAbs().maxCoeff();
            cout <<   "p = " << p << setw(15)
                 << "Max = "
                 << scientific << setprecision(3)
                 << err_max << endl;
        }
    }
    /* SAM_LISTING_END_1 */

    /* SAM_LISTING_BEGIN_4 */
    {
        auto u = [](double t) { return 2./M_PI*sqrt(t); };
        auto y = [](double t) { return t; };

        cout << "\n\nConvolution Quadrature, Implicit Euler\n"  << endl;
        for(int N=16; N<=2048; N*=2) {

            VectorXd grid = VectorXd::LinSpaced(N+1,0.,1.);
            VectorXd u_ex(N+1);
            for(int i=0; i<N+1; ++i) {
                u_ex(i) = u(grid(i));
            }
            
            VectorXd u_app = cq_ieul_abel(y, N);
            VectorXd diff  = u_ex - u_app;
            double err_max = diff.cwiseAbs().maxCoeff();
            cout <<   "N = " << N << setw(15)
                 << "Max = "
                 << scientific << setprecision(3)
                 << err_max << endl;
        }
        
        cout << "\n\nConvolution Quadrature, BDF-2\n"  << endl;
        for(int N=16; N<=2048; N*=2) {
            
            VectorXd grid = VectorXd::LinSpaced(N+1,0.,1.);
            VectorXd u_ex(N+1);
            for(int i=0; i<N+1; ++i) {
                u_ex(i) = 2./M_PI*sqrt(grid(i));
            }
            
            VectorXd u_app = cq_bdf2_abel(y, N);
            VectorXd diff  = u_ex - u_app;
            double err_max = diff.cwiseAbs().maxCoeff();
            cout <<   "N = " << N << setw(15)
                 << "Max = "
                 << scientific << setprecision(3)
                 << err_max << endl;
        }
    }
    /* SAM_LISTING_END_4 */

}
