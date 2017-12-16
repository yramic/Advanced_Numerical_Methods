#include <Eigen/Dense>
#include <unsupported/Eigen/FFT>
#include <cmath>
#include <iomanip>
#include <iostream>
#include "gauleg.hpp"

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
#if SOLUTION
    MatrixXd A = MatrixXd::Zero(p,p);
    VectorXd b = VectorXd::Zero(p);

    Eigen::RowVectorXd gauss_pts_p, gauss_wht_p;
    std::tie(gauss_pts_p,gauss_wht_p) = gauleg(0., 1., p);

    for(int i=1; i<=p; ++i) {

        for(int j=1; j<=p; ++j) {

            A(i-1,j-1) = 2.*sqrt(M_PI)*tgamma(1+j) / ((3.+2.*i+2.*j)*tgamma(3./2.+j)); // tgamma(1.+j) == j! is j is integer
        }

        for(int k=0; k<p; ++k) {

            double tk = gauss_pts_p(k);
            double wk = gauss_wht_p(k);

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

    return u;
#else // TEMPLATE
    // TODO: Find the unknown function u in the Abel integral equation with Galerkin discretization
#endif // TEMPLATE
}
/* SAM_LISTING_END_0 */


MatrixXd toeplitz_triangular(const VectorXd& c)
{
    size_t n = c.size();
    MatrixXd T = MatrixXd::Zero(n, n);
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
#if SOLUTION
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
#else // TEMPLATE
    // TODO: Find the unknown function u in the Abel integral equation with convolution quadrature (implicit Euler)
#endif // TEMPLATE
}
/* SAM_LISTING_END_2 */


VectorXcd pconvfft(const VectorXcd& u, const VectorXcd& x)
{
    FFT<double> fft;
    VectorXcd tmp = ( fft.fwd(u) ).cwiseProduct( fft.fwd(x) );
    return fft.inv(tmp);
}


VectorXcd myconv(const VectorXcd& h, const VectorXcd& x) {
  const long n = h.size();
  // Zero padding, cf. \eqref{eq:zeropad}
  VectorXcd hp(2*n - 1), xp(2*n - 1);
  hp << h, VectorXcd::Zero(n - 1);
  xp << x, VectorXcd::Zero(n - 1);
  // Periodic discrete convolution of length \Blue{$2n-1$}, \cref{cpp:pconffft}
  return pconvfft(hp, xp);
}


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
#if SOLUTION
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
#else // TEMPLATE
    // TODO: Find the unknown function u in the Abel integral equation with convolution quadrature (BDF-2)
#endif // TEMPLATE
}
/* SAM_LISTING_END_3 */


int main() {
    /* SAM_LISTING_BEGIN_1 */
#if SOLUTION
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
#else // TEMPLATE
    // TODO: Tabulate the max error of the Galerkin approximation scheme
#endif // TEMPLATE
    /* SAM_LISTING_END_1 */

    /* SAM_LISTING_BEGIN_4 */
#if SOLUTION
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
            double err_max = diff.cwiseAbs().maxCoeff();
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
            double err_max = diff.cwiseAbs().maxCoeff();
            cout <<   "N = " << N << setw(15)
                 << "Max = "
                 << scientific << setprecision(3)
                 << err_max << endl;
        }
    }
#else // TEMPLATE
    // TODO: Tabulate the max error of the convolution quadratures
#endif // TEMPLATE
    /* SAM_LISTING_END_4 */
}
