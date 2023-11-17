#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <Eigen/SparseLU>
#include <cmath>
#include <iomanip>
#include <iostream>

#include "absorbingboundarycondition.h"
using namespace Eigen;
using namespace std;


int main() {
  /* SAM_LISTING_BEGIN_3 */
#if SOLUTION
  auto g = [](double t) { return sin(M_PI * t); };
  // compute reference solution
  int M_ref = 4096;
  int N_ref = 4096;
  double h_ref = 1. / N_ref;
  double T = 1.0;
  VectorXd u_ref = AbsorbingBoundaryCondition::solve_IBVP(g, M_ref, N_ref, T);

  // compute H1 norm of reference solution
  double norm_u_ref = 0.;
  for (int i = 1; i <= N_ref; ++i) {
    norm_u_ref += pow((u_ref(i) - u_ref(i - 1)), 2);
  }
  norm_u_ref = sqrt(norm_u_ref / h_ref);

  cout << "\nConvergence wrt spatial discretisation" << endl;
  cout << "N\tRelative H1-error" << endl;
  for (int N = 16; N <= N_ref / 4; N *= 2) {
    double h = 1. / N;
    VectorXd u_tmp = AbsorbingBoundaryCondition::solve_IBVP(g, M_ref, N, T);
    double error = 0.;
    double ratio = N_ref / N;
    for (int i = 1; i <= N_ref; ++i) {
      int j = ceil(i / ratio);
      error += pow(
          (u_ref(i) - u_ref(i - 1)) / h_ref - (u_tmp(j) - u_tmp(j - 1)) / h, 2);
    }
    error = sqrt(error * h_ref) / norm_u_ref;
    cout << N << "\t" << scientific << setprecision(10) << error << endl;
  }

  cout << "\nConvergence wrt time discretisation" << endl;
  cout << "M\tRelative H1-error" << endl;
  for (int M = 16; M <= M_ref / 4; M *= 2) {
    VectorXd u_tmp = AbsorbingBoundaryCondition::solve_IBVP(g, M, N_ref, T);
    double error = 0.;
    for (int i = 1; i <= N_ref; ++i) {
      error += pow((u_ref(i) - u_ref(i - 1)) - (u_tmp(i) - u_tmp(i - 1)), 2);
    }
    error = sqrt(error / h_ref) / norm_u_ref;
    cout << M << "\t" << scientific << setprecision(10) << error << endl;
  }

  // Call python script
  std::system("python3 " CURRENT_SOURCE_DIR "/plot.py " CURRENT_BINARY_DIR
              "/u_ref.txt " CURRENT_BINARY_DIR "/u.png");
#else   // TEMPLATE
  // TODO: Tabulate the H1-error of the Galerkin discretization + convolution
  // quadrature
#endif  // TEMPLATE
  /* SAM_LISTING_END_3 */
}
