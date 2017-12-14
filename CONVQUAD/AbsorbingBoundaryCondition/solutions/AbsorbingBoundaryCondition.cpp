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
/* SAM_LISTING_BEGIN_0 */
template<typename FUNC>
VectorXd poly_spec_abel(const FUNC& y, size_t p, double tau)
{
    MatrixXd A = MatrixXd::Zero(p,p);
    VectorXd b = VectorXd::Zero(p);

    const double* gauss_pts_p = getGaussPoints(p);
    const double* gauss_wht_p = getGaussWeights(p);

    for(int i=1; i<=p; ++i) {

        for(int j=1; j<=p; ++j) {

            A(i-1,j-1) = 2.*sqrt(M_PI)*tgamma(1.+j) / ((3.+2.*i+2.*j)*tgamma(3./2.+j)); // tgamma(1.+j) == j! is j is integer
        }

        for(int k=0; k<p; ++k) {

            double tk = 0.5 * (gauss_pts_p[k] + 1.);
            double wk = 0.5 *  gauss_wht_p[k];

            b(i-1) += wk * pow(tk,i) * y(tk);
        }
    }

    VectorXd x = A.colPivHouseholderQr().solve(b);

    size_t N = round(1./tau);
    VectorXd grid = VectorXd::LinSpaced(N+1, 0., 1.);
    VectorXd u    = VectorXd::Zero(N+1);

    for(int i=0; i<N+1; ++i) {
        for(int j=1; j<=p; ++j) {
            u(i) += x(j-1) * pow(grid(i),j);
        }
    }

    return grid;
}
/* SAM_LISTING_END_0 */


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
/* SAM_LISTING_BEGIN_2 */
template<typename FUNC>
VectorXd cq_ieul_abel(const FUNC& y, size_t N)
{
    VectorXcd w = VectorXd::Zero(N+1);
    double tau = 1./N;

    int p = 8; // order of quadrature rule
    const double* gauss_pts_p = getGaussPoints(p);
    const double* gauss_wht_p = getGaussWeights(p);

    for(int i=0; i<p; ++i) {

        // integrate on semi-circumference centered in (1,0) with unitary radius:
        complex<double> ti = 1. + exp( complex<double>(0.,-0.5*M_PI*gauss_pts_p[i]) ); // M\_PI/2 - 0.5*(gauss\_pts\_p[i]+1.)*M\_PI
                     double  wi = 0.5 * M_PI * gauss_wht_p[i]; // change of integration domain to semi-circumference

        for(int j=0; j<N+1; ++j) {
            w(j) += wi / (sqrt(ti)*pow(1.-tau*ti,j+1));
        }

        // integrate on segment from (+1,-1) to (+1,+1):
        ti = complex<double>(1.,gauss_pts_p[i]);
        wi = gauss_wht_p[i];

        for(int j=0; j<N+1; ++j) {
            w(j) += wi / (sqrt(ti)*pow(1.-tau*ti,j+1));
        }
    }

    w *= tgamma(0.5)*tau / complex<double>(0.,2.*M_PI);

    // Solve the convolution quadrature:

    VectorXd  grid = VectorXd::LinSpaced(N+1,0.,1.);
    VectorXcd y_N(N+1);
    for(int i=0; i<N+1; ++i) {
        y_N(i) = complex<double>(y(grid(i)),0.);
    }

    MatrixXcd T = toeplitz_triangular(w);
    VectorXcd u = T.triangularView<Lower>().solve(y_N);

    return u.real();
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
    VectorXcd w = VectorXd::Zero(N+1);
    double tau = 1./N;

    int p = 8; // order of quadrature rule
    const double* gauss_pts_p = getGaussPoints(p);
    const double* gauss_wht_p = getGaussWeights(p);

    for(int i=0; i<p; ++i) {

        // integrate on semi-circumference centered in (1,0) with unitary radius:
        complex<double> ti = 1. + exp( complex<double>(0.,-0.5*M_PI*gauss_pts_p[i]) ); // M\_PI/2 - 0.5*(gauss\_pts\_p[i]+1.)*M\_PI
                     double  wi = 0.5 * M_PI * gauss_wht_p[i]; // change of integration domain to semi-circumference

        for(int j=0; j<N+1; ++j) {
            w(j) += wi / (sqrt(ti)*sqrt(1.+tau*ti)) * (1./pow(2.-sqrt(1.+2.*tau*ti),j+1) - 1./pow(2.+sqrt(1.+2.*tau*ti),j+1));
        }

        // integrate on segment from (+1,-1) to (+1,+1):
        ti = complex<double>(1.,gauss_pts_p[i]);
        wi = gauss_wht_p[i];

        for(int j=0; j<N+1; ++j) {
            w(j) += wi / (sqrt(ti)*sqrt(1.+tau*ti)) * (1./pow(2.-sqrt(1.+2.*tau*ti),j+1) - 1./pow(2.+sqrt(1.+2.*tau*ti),j+1));
        }
    }

    w *= tgamma(0.5)*tau / complex<double>(0.,2.*M_PI);

    // Solve the convolution quadrature:

    VectorXd  grid = VectorXd::LinSpaced(N+1,0.,1.);
    VectorXcd y_N(N+1);
    for(int i=0; i<N+1; ++i) {
        y_N(i) = complex<double>(y(grid(i)),0.);
    }

    MatrixXcd T = toeplitz_triangular(w);
    VectorXcd u = T.triangularView<Lower>().solve(y_N);

    return u.real();
}
/* SAM_LISTING_END_3 */


int main() {
    /* SAM_LISTING_BEGIN_1 */
    {
        auto y = [](double t) { return t; };

        double tau = 0.01;
        size_t N = round(1./tau);
        VectorXd grid = VectorXd::LinSpaced(N+1,0.,1.);
        VectorXd u_ex(N+1);
        for(int i=0; i<N+1; ++i) {
            u_ex(i) = 2./M_PI*sqrt(grid(i));
        }

        cout << "Problem 3.1.g" << endl;
        for(int p=2; p<=32; p*=2) {
            VectorXd u_app = poly_spec_abel(y, p, tau);
            VectorXd diff  = u_ex - u_app;
            double err_max  = diff.cwiseAbs().maxCoeff();
            cout <<   "p = " << p << setw(15)
                      << "Max = "
                      << scientific << setprecision(3)
                      << err_max << endl;
        }
    }
    /* SAM_LISTING_END_1 */

    /* SAM_LISTING_BEGIN_4 */
    {
        auto y = [](double t) { return t; };

        cout << "Problem 3.1.l"  << endl;
        cout << "Implicit Euler" << endl;
        for(int N=10; N<=1280; N*=2) {

            VectorXd grid = VectorXd::LinSpaced(N+1,0.,1.);
            VectorXd u_ex(N+1);
            for(int i=0; i<N+1; ++i) {
                u_ex(i) = 2./M_PI*sqrt(grid(i));
            }

            VectorXd u_app = cq_ieul_abel(y, N);
            VectorXd diff  = u_ex - u_app;
            double err_max  = diff.cwiseAbs().maxCoeff();
            cout <<   "N = " << N << setw(15)
                      << "Max = "
                      << scientific << setprecision(3)
                      << err_max << endl;
        }

        cout << "BDF-2" << endl;
        for(int N=10; N<=1280; N*=2) {

            VectorXd grid = VectorXd::LinSpaced(N+1,0.,1.);
            VectorXd u_ex(N+1);
            for(int i=0; i<N+1; ++i) {
                u_ex(i) = 2./M_PI*sqrt(grid(i));
            }

            VectorXd u_app = cq_bdf2_abel(y, N);
            VectorXd diff  = u_ex - u_app;
            double err_max  = diff.cwiseAbs().maxCoeff();
            cout <<   "N = " << N << setw(15)
                      << "Max = "
                      << scientific << setprecision(3)
                      << err_max << endl;
        }
    }
    /* SAM_LISTING_END_4 */
}
