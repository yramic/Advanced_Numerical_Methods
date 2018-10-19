#include <iostream>
#include <fstream>
#include <istream>
#include <cmath>
#include <Eigen/Dense>
#include <boost/math/special_functions/bessel.hpp> // For Bessel Function
#include "chebyshev_gauss_quadrature.hpp"
#include <boost/math/special_functions/chebyshev.hpp>

using namespace Eigen;

#if SOLUTION
#include "PeriodicTrapezoidalQR.hpp"
#endif // SOLUTION

/* @brief Compute matrix M using analytic expression.
 *
 * \param[in] N Discretization parameter indicating number of basis functions.
 */
 /* SAM_LISTING_BEGIN_0 */
MatrixXd computeM (unsigned int N) {
  Eigen::MatrixXd M(N,N);
  M.setZero();

#if SOLUTION
// Matrix assembly
for (unsigned int i = 0 ; i < N ; ++i)
  M(i,i) = M_PI*M_PI/2./(i+1);

#else // TEMPLATE
  // TODO: Compute M
#endif // TEMPLATE
  return M;
}
/* SAM_LISTING_END_0 */

/* @brief Compute right hand side using periodic trapezoidal rule (2N points).
 *
 * \tparam FUNC Type for the Dirichlet B.C. supporting : g(const double&)
 * \param[in] g Dirichlet boundary condition; Signature : double(const double&)
 * \param[in] N Discretization parameter indicating number of basis functions.
 */
 /* SAM_LISTING_BEGIN_1 */
template <typename FUNC>
VectorXd computeG(const FUNC& g, int N){
  // Initialize right hand side vector
  VectorXd RHS(N);  RHS.setZero();
  #if SOLUTION
  // Get Chebyshev Gauss Quadrature points and weight
  unsigned int order = 100;
  Eigen::RowVectorXd points;
  double weight;
  std::tie(points,weight) = ChebyshevGaussQuad(order);
  // Fill vector entries
  for (unsigned int i = 0 ; i < N ; ++i) {
    // Evaluate integral by using the quadrature
    for (unsigned int qp = 0 ; qp < order ; ++qp) {
      double x = points(qp);
      RHS(i) += weight * g(x)*boost::math::chebyshev_t(i+1, x);
    }// end iteration over quadrature points
  }

  #else // TEMPLATE
    // TODO: Compute RHS
  #endif // TEMPLATE
  return RHS;
}
/* SAM_LISTING_END_1 */

/* @brief Build and solve boundary integral equation V rho = g
 *
 * \tparam FUNC Type for the Dirichlet B.C. supporting : g(const double&)
 * \param[in] g Dirichlet boundary condition; Signature : double(const double&)
 * \param[in] N Discretization parameter indicating number of basis functions.
 */
 /* SAM_LISTING_BEGIN_2 */
template <typename FUNC>
VectorXd solveBIE(const FUNC& g, int N){
#if SOLUTION
  // Assemble Galerkin Matrix M
  Eigen::MatrixXd M = computeM(N);
  // Build RHS
  VectorXd RHS = computeG(g, N);
  // Use direct solver
  MatrixXd LHS = M.eval();
  VectorXd sol = LHS.lu().solve(RHS);
#else // TEMPLATE
    // TODO: Build BIE system and solve it
#endif // TEMPLATE

  return sol;
}
/* SAM_LISTING_END_2 */

/* @brief Reconstruct function UN from its coefficients and evaluate it at t.
 *
 * \param[in] coeffs coefficients of UN
 * \param[in] t evaluation point [-1,1]
 */
double reconstructRho(const VectorXd& coeffs, double t){
  int N = coeffs.rows();
  double rho_N = 0.;
  for (unsigned int i = 0 ; i < N ; ++i)
    rho_N += coeffs(i) * boost::math::chebyshev_t(i+1, t);
  return rho_N/std::sqrt(1-t*t);
}


/* @brief Compute L2 norm of UN from its coefficients using periodic trapezoidal
 *        rule (2N points).
 *
 * \param[in] coeffs coefficients of UN
 */
 /*
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
}*/

/* @brief Compute L2 norm of a function on the domain [-1,1] using Gauss Quadrature
 *
 * \tparam Function Template type for the function whose L2 norm is to be
 *                  evaluated. Must support evaluation : double Function(double)
 * \param[in] function The function for which the L2 norm has to be evaluated
 * \param[in] N Number of Gauss Quadrature points to be used for norm evaluation
 * \param[in] coeffs coefficients of UN
 */
template<typename Function>
double l2norm(const Function& function,const int& N) {
  double norm = 0.;
  // Get quadrature points and weight for Gauss Legendre Quadrature
  unsigned int order = 200;
  Eigen::RowVectorXd points;
  double weight;
  std::tie(points,weight) = ChebyshevGaussQuad(order);
  for(int qp=0; qp<order; qp++){
    auto z = points(qp);
    // evaluate
    double rho = function(z);
    norm += weight*rho*rho;
  }
  return std::sqrt(norm);
}

double l2wnorm(const VectorXd& coeffs) {
  int N = coeffs.rows();

  // Lambda function to evaluate rho_N * w
  std::function<double(const double&)> rho_N_w = [&](const double& x) {
    double rho_N = 0.;
    for (unsigned int i = 0 ; i < N ; ++i)
      rho_N += coeffs(i) * boost::math::chebyshev_t(i+1, x);
    return rho_N;
   };

   // Lambda function to evaluate rho1_exact * w
   std::function<double(const double&)> rho1_exact_w = [&](const double& x) {
     unsigned int trunc_num = 100;
     double rho_exact = 0.;
     for (unsigned int n = 0 ; n < trunc_num ; ++n)
      rho_exact += 4*n/M_PI*boost::math::cyl_bessel_j(n,2*M_PI)*sin(n*M_PI/2.)*boost::math::chebyshev_t(2*n, x);
     return rho_exact;
    };
  double norm = 0.;
  // Get quadrature points and weight for Gauss Chebyshev Quadrature
  unsigned int order = 100;
  Eigen::RowVectorXd points;
  double weight;
  std::tie(points,weight) = ChebyshevGaussQuad(order);
  // Integration using quadrature
  for(int qp=0; qp<order; qp++){
    auto x = points(qp);
    double integrand = rho_N_w(x)-rho1_exact_w(x); // error
    norm += weight*integrand*integrand; // calculating norm of the error
  }
  return std::sqrt(norm);
}

double rho1_exact(double x) {
  //Truncating the infinite sum for exact solution of rho1(x)
  unsigned int trunc_num = 1000;
  double rho_exact = 0.;
  for (unsigned int n = 0 ; n < trunc_num ; ++n)
   rho_exact += 4*n/M_PI*boost::math::cyl_bessel_j(n,2*M_PI)*sin(n*M_PI/2.)*boost::math::chebyshev_t(2*n, x);
  return rho_exact/std::sqrt(1-x*x);
}

int main() {
  std::cout << "Test for source term g1(x) = sin(2*Pi*x)" << std::endl;
  std::cout << "N" << std::setw(10)  << "L2error" << std::setw(12) << "L2Werror" << std::endl;
  std::cout << "############################" << std::endl;
  double test_pt = 0.33;

  std::function<double(const double&)> g1 = [](const double& x) {
     return sin(2*M_PI*x);
   };
  // Varying the discretization parameter
  for (unsigned int N = 2 ; N < 50 ; ++N) {
    // Solving the BIE to get the coefficients
    //unsigned int N = 8;
    Eigen::MatrixXd M = computeM(N);
    Eigen::VectorXd coeffs = solveBIE(g1,N);
    //std::cout << "Evaluated coefficients: " << coeffs << std::endl;
    // Function to evaluate reconstructed rho at the point x
    std::function<double(const double&)> rho1_N = [&](const double& x) {
       return reconstructRho(coeffs,x);
     };
    // Function to evaluate rho_exact-reconstructed rho at the point x
    std::function<double(const double&)> delta = [&](const double& x) {
       return rho1_exact(x)-reconstructRho(coeffs,x);
     };

    // L2 norm error
    double l2error = l2norm(delta,N);
    double l2werror = l2wnorm(coeffs);
    std::cout << N << std::setw(10) << l2error << std::setw(12) << l2werror << std::endl;

  }
  return 0;
}
