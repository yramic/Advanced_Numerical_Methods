#include <Eigen/Dense>
#include <cmath>
#include <iomanip>
#include <iostream>

using namespace Eigen;
using namespace std;


/* @brief Find the unknown function u in the Abel integral equation
 * using convolution quadrature (BDF-2)
 * \param y Template function for the right-hand side
 * \param N Number of discretization steps
 * \\return Values of u from convolution quadrature
 */
/* SAM_LISTING_BEGIN_0 */
template<typename FUNC>
VectorXd cq_bdf2_abel(const FUNC& y, size_t N)
{
//#if SOLUTION
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
#else // TEMPLATE
    // TODO: Find the unknown function u in the Abel integral equation with convolution quadrature (BDF-2)
#endif // TEMPLATE
}
/* SAM_LISTING_END_0 */


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
#else // TEMPLATE
    // TODO: Tabulate the max error of the convolution quadratures
#endif // TEMPLATE
    /* SAM_LISTING_END_4 */
}
