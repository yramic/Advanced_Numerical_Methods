#include <iostream>
#include <fstream>
#include <istream>
#include <cmath>
#include <Eigen/Dense>
#include <boost/math/special_functions/bessel.hpp> // For Bessel Function
#include "gauleg.hpp"

using namespace Eigen;

#if SOLUTION
#include "PeriodicTrapezoidalQR.hpp"
#endif // SOLUTION

/* @brief Compute matrix M using analytic expression.
 *
 * \param[in] N Discretization parameter indicating number of basis functions.
 */
MatrixXd computeM (unsigned int N) {
  Eigen::MatrixXd M(N,N);
  M.setZero();

#if SOLUTION
// Matrix assembly
for (unsigned int i = 0 ; i < N ; ++i)
  M(i,i) = M_PI/4./(i+1);

#else // TEMPLATE
  // TODO: Compute M
#endif // TEMPLATE
  return M;
}

/* @brief Compute right hand side using periodic trapezoidal rule (2N points).
 *
 * \tparam FUNC Type for the Dirichlet B.C. supporting : g(const double&)
 * \param[in] g Dirichlet boundary condition; Signature : double(const double&)
 * \param[in] N Discretization parameter indicating number of basis functions.
 */
template <typename FUNC>
VectorXd computeG(const FUNC& g, int N){
  // Initialize right hand side vector
  VectorXd RHS(N);  RHS.setZero();
  #if SOLUTION
  // Get periodic Trapezoidal Quadrature points and weight
  //VectorXd TR_points(2*N);  double TR_w;
  //std::tie(TR_points,TR_w) = PeriodicTrapRule(2*N);
  Eigen::RowVectorXd weights,points;
  std::tie(points,weights) = gauleg(0,M_PI,2*N);
  // Fill vector entries
  for (unsigned int i = 0 ; i < N ; ++i) {
    // Evaluate integral by using the quadrature
    for (unsigned int qp = 0 ; qp < 2*N ; ++qp) {
      double theta = points(qp);
      RHS(i) += weights(qp) * g(cos(theta))*cos((i+1)*theta);
    }// end iteration over quadrature points
  }

  #else // TEMPLATE
    // TODO: Compute RHS
  #endif // TEMPLATE
  return RHS;
}

/* @brief Build and solve boundary integral equation V rho = g
 *
 * \tparam FUNC Type for the Dirichlet B.C. supporting : g(const double&)
 * \param[in] g Dirichlet boundary condition; Signature : double(const double&)
 * \param[in] N Discretization parameter indicating number of basis functions.
 */
template <typename FUNC>
VectorXd solveBIE(const FUNC& g, int N){
#if SOLUTION
  std::cout << " Assemble M " << std::endl;
  Eigen::MatrixXd M = computeM(N);
  // Build RHS
  std::cout << " Assemble RHS" << std::endl;
  VectorXd RHS = computeG(g, N);
  // Use direct solver
  std::cout << " Solve " << std::endl;
  MatrixXd LHS = M.eval();
  VectorXd sol = LHS.lu().solve(RHS);
  std::cout << " Done " << std::endl;
#else // TEMPLATE
    // TODO: Build BIE system and solve it
#endif // TEMPLATE

  return sol;
}

/* @brief Reconstruct function UN from its coefficients and evaluate it at t.
 *
 * \param[in] coeffs coefficients of UN
 * \param[in] t evaluation point [-1,1]
 */
double reconstructRho(const VectorXd& coeffs, double t){
  int N = coeffs.rows();
  double rho_N = 0.;
  for (unsigned int i = 0 ; i < N ; ++i) {
    rho_N += coeffs(i) * cos((i+1) * acos(t));
  }
  return rho_N/std::sqrt(1-t*t);
}

/* @brief Compute L2 norm of UN from its coefficients using periodic trapezoidal
 *        rule (2N points).
 *
 * \param[in] coeffs coefficients of UN
 */
double L2norm(const VectorXd& coeffs){
  double norm = 0.;
#if SOLUTION
  int N = coeffs.rows();
  // Get quadrature points and weight
  VectorXd TR_points(2*N);  double TR_w;
  std::tie(TR_points,TR_w) = PeriodicTrapRule(2*N);
  for(int qp=0; qp<2*N; qp++){
    auto z = TR_points(qp);
    // evaluate
    double rho = reconstructRho(coeffs, z);
    norm += TR_w*rho*rho;
  }
  #else // TEMPLATE
    // TODO: reconstruct the function and compute its L2norm
  #endif // TEMPLATE

  return std::sqrt(norm);
}

int main() {
  std::cout << "Test for source term g1(x) = sin(2*Pi*x)" << std::endl;
  double test_pt = 0.33;
   //Truncating the infinite sum for exact solution of rho1(x)
   unsigned int trunc_num = 5000;
   double rho1_exact = 0.;
  for (unsigned int n = 0 ; n < trunc_num ; ++n) {
    rho1_exact += 4*n/M_PI*boost::math::cyl_bessel_j(n,2*M_PI)*sin(n*M_PI/2.)*cos(2*n*acos(test_pt));
  }
  std::function<double(const double&)> g1 = [](const double& x) {
     return sin(2*M_PI*x);
   };
  // Varying the discretization parameter
  for (unsigned int N = 2 ; N < 20 ; ++N) {
  //unsigned int N = 10;
    /*Eigen::VectorXd trunc_coeffs(N);
    trunc_coeffs.setZero();
    for (unsigned int n = 0 ; n < N ; ++n) {
      if ((n+1)%2==1)
        trunc_coeffs(n) = 0;
      else
        trunc_coeffs(n) = 2.*(n+1)/M_PI*boost::math::cyl_bessel_j((n+1)/2.,2*M_PI)*sin((n+1)/2.*M_PI/2.);
    }*/
    Eigen::MatrixXd M = computeM(N);
    //std::cout << "M is \n" << M << std::endl;
    Eigen::VectorXd coeffs = solveBIE(g1,N);
    /*std::cout << "Solved coefficients" << std::endl;
    std::cout << coeffs << std::endl;
    std::cout << "Computed coefficients" << std::endl;
    std::cout << trunc_coeffs << std::endl;
    double error = L2norm(coeffs-trunc_coeffs);
    std::cout << "N = " << N << " error = " << error << std::endl;*/
    double rho1_N = reconstructRho(coeffs,test_pt);
    //std::cout << "Exact = " << rho1_exact << std::endl;
    //std::cout << "Computed = " << rho1_N << std::endl;
    double error = fabs(rho1_exact-rho1_N);
    std::cout << "N = " << N << " error = " << error << std::endl;
  }

  return 0;
}
