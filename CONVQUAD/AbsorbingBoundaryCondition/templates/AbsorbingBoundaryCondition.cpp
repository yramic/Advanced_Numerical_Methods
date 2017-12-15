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
    // TODO: Find the unknown function u at final time t = 1 in the evolution problem
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
