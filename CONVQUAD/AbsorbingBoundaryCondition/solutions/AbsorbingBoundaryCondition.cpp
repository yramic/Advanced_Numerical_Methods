#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <Eigen/SparseCholesky>
#include <cmath>
#include <iomanip>
#include <iostream>

using namespace Eigen;
using namespace std;


VectorXd pconv(const VectorXd& u, const VectorXd& x) {
  using idx_t = VectorXd::Index; // may be unsigned !
  const idx_t n = x.size();
  VectorXd z = VectorXd::Zero(n);
  // Need signed indices when differences are formed
  for (long k = 0; k < n; ++k) {
      for (long j = 0; j < n; ++j) {
          long ind = (k - j < 0 ? n + k - j : k - j);
          z(k) += u(ind)*x(j);
      }
  }
  return z;
}


/* @brief Find the unknown function u at final time t = 1 in the evolution problem
 * using Galerkin discretization and convolution quadrature (BDF-2)
 * \param phi Template function for the right-hand side
 * \param M Number of discretization steps in time
 * \param N Number of discretization steps in space
 * \param p Order of quadrature rule
 * \\return Values of u at final time t = 1
 */
/* SAM_LISTING_BEGIN_0 */
template<typename FUNC>
VectorXd solveABC(const FUNC& phi, size_t M, size_t N, int p)
{
    double h = 1./(N-1);
    VectorXd grid = VectorXd::LinSpaced(N,0.,1.);

    SparseMatrix<double> A(N,N);
    std::vector<Eigen::Triplet<double> > triplets(3*N-2);
    triplets[0] = Eigen::Triplet<double>(0,0, 1./h - 1./(h*h*pow(M_PI,3))*(2. + (M_PI*M_PI*h*h-2.)*cos(M_PI*h) - 2.*M_PI*h*sin(M_PI*h)) );
    triplets[1] = Eigen::Triplet<double>(1,0,-1./h + 1./(h*h*pow(M_PI,3))*(2. - 2.*cos(M_PI*h) - M_PI*h*sin(M_PI*h)) );
    for(int i=1; i<N-1; ++i) {
        triplets[2+3*i  ] = Eigen::Triplet<double>(i-1,i,-1./h + 1./(h*h*pow(M_PI,3))*(2.*(cos(M_PI*grid(i-1)) - cos(M_PI*grid(i))) - M_PI*h*(sin(M_PI*grid(i-1)) + sin(M_PI*grid(i)))) );
        triplets[2+3*i+1] = Eigen::Triplet<double>(i,  i, 2./h + 1./(h*h*pow(M_PI,3))*((M_PI*M_PI*h*h-4.)*(cos(M_PI*grid(i-1)) - cos(M_PI*grid(i+1))) + 4.*M_PI*h*(sin(M_PI*grid(i-1)) + sin(M_PI*grid(i+1)))) );
        triplets[2+3*i+2] = Eigen::Triplet<double>(i+1,i,-1./h + 1./(h*h*pow(M_PI,3))*(2.*(cos(M_PI*grid(i)) - cos(M_PI*grid(i+1))) - M_PI*h*(sin(M_PI*grid(i)) + sin(M_PI*grid(i+1)))) );
    }
    triplets[3*N-4] = Eigen::Triplet<double>(N-2,N-1,-1./h + 1./(h*h*pow(M_PI,3))*(2.*cos(M_PI*(1.-h)) - 2. - M_PI*h*sin(M_PI*(1.-h))) );
    triplets[3*N-3] = Eigen::Triplet<double>(N-1,N-1, 1./h + 1./(h*h*pow(M_PI,3))*((M_PI*M_PI*h*h-2.)*cos(M_PI*(1.-h)) - 2. + 2.*M_PI*h*sin(M_PI*(1.-h))) );

    A.setFromTriplets(triplets.begin(), triplets.end());

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
    Aw.coeffRef(N-1,N-1) += w(0);
    SimplicialLLT<SparseMatrix<complex<double> > > solver;
    solver.compute(Aw);

    MatrixXcd u = VectorXcd::Zero(N,M+1); // first column is at time t = 0
    for(int i=1; i<M+1; ++i) {

        complex<double> rhs_N = 0.;
        for(int l=0; l<i; ++l) {
            rhs_N += w(i-l) * u(N-1,l);
        }

        VectorXcd phi_i(N);
        for(int j=0; j<N; ++j) {
            phi_i(j) = (complex<double>)phi(grid(j),i*tau);
        }

        VectorXcd rhs = phi_i; rhs(N-1) -= rhs_N;

        u.col(i) = solver.solve(rhs);
    }

    return u.col(M).real();
}
/* SAM_LISTING_END_0 */


int main() {
//    /* SAM_LISTING_BEGIN_1 */
//#if SOLUTION
//    {
//        auto y = [](double t) { return t; };

//        double tau = 0.01;
//        size_t N = round(1./tau);
//        VectorXd grid = VectorXd::LinSpaced(N+1,0.,1.);
//        VectorXd u_ex(N+1);
//        for(int i=0; i<N+1; ++i) {
//            u_ex(i) = 2./M_PI*sqrt(grid(i));
//        }

//        cout << "Problem 3.1.g" << endl;
//        for(int p=2; p<=32; p*=2) {
//            VectorXd u_app = poly_spec_abel(y, p, tau);
//            VectorXd diff  = u_ex - u_app;
//            double err_max  = diff.cwiseAbs().maxCoeff();
//            cout <<   "p = " << p << setw(15)
//                      << "Max = "
//                      << scientific << setprecision(3)
//                      << err_max << endl;
//        }
//    }
//#else // TEMPLATE
//    // TODO: Tabulate the max error of the Galerkin approximation scheme
//#endif // TEMPLATE
//    /* SAM_LISTING_END_1 */

//    /* SAM_LISTING_BEGIN_4 */
//#if SOLUTION
//    {
//        auto y = [](double t) { return t; };

//        cout << "Problem 3.1.l"  << endl;
//        cout << "Implicit Euler" << endl;
//        for(int N=10; N<=1280; N*=2) {

//            VectorXd grid = VectorXd::LinSpaced(N+1,0.,1.);
//            VectorXd u_ex(N+1);
//            for(int i=0; i<N+1; ++i) {
//                u_ex(i) = 2./M_PI*sqrt(grid(i));
//            }

//            VectorXd u_app = cq_ieul_abel(y, N);
//            VectorXd diff  = u_ex - u_app;
//            double err_max  = diff.cwiseAbs().maxCoeff();
//            cout <<   "N = " << N << setw(15)
//                      << "Max = "
//                      << scientific << setprecision(3)
//                      << err_max << endl;
//        }

//        cout << "BDF-2" << endl;
//        for(int N=10; N<=1280; N*=2) {

//            VectorXd grid = VectorXd::LinSpaced(N+1,0.,1.);
//            VectorXd u_ex(N+1);
//            for(int i=0; i<N+1; ++i) {
//                u_ex(i) = 2./M_PI*sqrt(grid(i));
//            }

//            VectorXd u_app = cq_bdf2_abel(y, N);
//            VectorXd diff  = u_ex - u_app;
//            double err_max  = diff.cwiseAbs().maxCoeff();
//            cout <<   "N = " << N << setw(15)
//                      << "Max = "
//                      << scientific << setprecision(3)
//                      << err_max << endl;
//        }
//    }
//#else // TEMPLATE
//    // TODO: Tabulate the max error of the convolution quadratures
//#endif // TEMPLATE
//    /* SAM_LISTING_END_4 */
}
