//// 
//// Copyright (C) 2017 SAM (D-MATH) @ ETH Zurich
//// Author(s): curzuato < > 
//// Contributors:  dcasati 
//// This file is part of the AdvNumCSE repository.
////
#include <iostream>
#include <fstream>
#include <istream>
#include <cmath>
#include <Eigen/Dense>
#include <iomanip>

using namespace Eigen;




/* @brief Compute matrix A-M using analytic expression.
 * \param[in] N Discretization parameter indicating number of basis functions.
 */
MatrixXd computeAminusM(int N){
  // A-M is a diagonal matrix of size 2N+1 X 2N+1
  DiagonalMatrix<double, Dynamic> AminusM(2*N+1);

  // TODO: Compute A-M
  return AminusM;
}




//----------------------------------------------------------------------------
/* @brief Compute matrix M using periodic trapezoidal rule (2N+2 points).
 * \tparam PARAM Type for gamma which supports evaluation operator: gamma(double)
 * \param[in] gamma Function that takes a double and returns a 2d vector
 *                  corresponding to the parametrized curve.
 * \param[in] N Discretization parameter indicating number of basis functions.
 */
template <typename PARAM>
MatrixXd computeM(const PARAM& gamma, int N){

  // Assemble big matrix
  MatrixXd M(2*N+1, 2*N+1);

    // TODO: Compute M using quadrature

  return M;
}


//----------------------------------------------------------------------------
/* @brief Compute right hand side using periodic trapezoidal rule (2N points).
 * \tparam PARAM Type for gamma which supports evaluation operator: gamma(double)
 * \tparam FUNC Type for function g which supports evaluation operator: g(double)
 * \param[in] gamma Function that takes a double and returns a 2d vector
 *                  corresponding to the parametrized curve.
 * \param[in] g Right hand side function (takes 2d points and returns a double).
 * \param[in] N Discretization parameter indicating number of basis functions.
 */
template <typename PARAM, typename FUNC>
VectorXd computeG(const PARAM& gamma, const FUNC& g, int N){
  // Initialize right hand side vector
  VectorXd RHS(2*N+1);  RHS.setZero();
    // TODO: Compute RHS
  return RHS;
}


//----------------------------------------------------------------------------
/* @brief Build and solve boundary integral equation V rho = g
 * \tparam PARAM Type for gamma which supports evaluation operator: gamma(double)
 * \tparam FUNC Type for function g which supports evaluation operator: g(double)
 * \param[in] gamma Function that takes a double and returns a 2d vector
 *                  corresponding to the parametrized curve.
 * \param[in] g Right hand side function (takes 2d points and returns a double).
 * \param[in] N Discretization parameter indicating number of basis functions.
 */
template <typename PARAM, typename FUNC>
VectorXd solveBIE(const PARAM& gamma, const FUNC& g, int N){
    // TODO: Build BIE system and solve it
  VectorXd sol(2*N+1);

  return sol;
}


//----------------------------------------------------------------------------


//----------------------------------------------------------------------------
/* @brief Compute L2 norm of UN from its coefficients using periodic trapezoidal
 *        rule (2N points).
 * \tparam PARAMDER Type for gammaprime which supports evaluation operator:
 *                  gammaprime(double)
 * \param[in] gammaprime Function that takes a double and returns a 2d vector
 *                       corresponding to the derivative of the curve's
 *                       parametrization.
 * \param[in] coeffs coefficients of UN
 */
template <typename PARAMDER>
double L2norm(const PARAMDER& gammaprime, const VectorXd& coeffs){
  double res=0.;
    // TODO: reconstruct the function and compute its L2norm

  return std::sqrt(res);
}





int main() {

  int N = 10;
  //----------------------------------------------------------------------------
  std::cout << "===========  Test integration errors ===========" << std::endl;
  double Qint1 = 0., Qintcos = 0., Qintlogcos2 = 0.;
  /* TODO: You may test your quadrature by integrating 1, cos(t) and
     log(cos(t)^2 + 1) */
  double intlogcos2ex = 2.*M_PI*log((sqrt(2)+1)/((sqrt(2)-1)*4));
  std::cout << " int 1                : " << fabs (Qint1 - 2*M_PI)
	    << "\n int cosx             : " << fabs(Qintcos)
	    << "\n int log(cosx*cosx+1) : " << fabs(Qintlogcos2 - intlogcos2ex)
	    << std::endl;
  std::cout << "================================================" << std::endl
	    << std::endl;


  //----------------------------------------------------------------------------
  std::cout << "============  Test Coefficients for S(t) = (cos(t), sin(t))  "
	    << "============="  << std::endl;
  N = 4;
  // Assigning a Lambda expression to the function S
  std::function<Vector2d(const double&)> S = [](const double& t){
    Vector2d res;
    res << cos(t) , sin(t);
    return res;
  };
  // Assigning a Lambda expression to the function Sprime
  std::function<Vector2d(const double&)> Sprime = [](const double& t){
    Vector2d res;
    res << -sin(t),  cos(t);
    return res;
  };


  // TODO: You may test your computation of the coefficients for gamma(t) and the
  // difference (gamma(0.1)-gamma(0))/||S_hat(0.1)-S_hat(0)||.
  Vector2d diffS; diffS.setZero();
  // Calculating gamma(s)-gamma(t) exactly
  Vector2d exDiffS; exDiffS<< cos(0.1) - cos(0), sin(0.1);
  std::cout << "(S(0.1)-S(0))/||S(0.1)-S(0)|| = " << std::endl;
  std::cout << "Calculated" << std::setw(20) << "Exact" <<  std::endl;
  std::cout << diffS.transpose()(0) << ", " << diffS.transpose()(1) << std::setw(16)
      // Calculating (gamma(s)-gamma(t))/||S_hat(s)-S_hat(t)|| ; gamma = S = S_hat
	    << exDiffS(0)/exDiffS.norm() << ", " << exDiffS(1)/exDiffS.norm()
	    << std::endl;
  std::cout << "=============================================================="
	    << "============" << std::endl << std::endl;


  //----------------------------------------------------------------------------
  std::cout << "====================  Test Coefficients for gamma(t)  "
	    << "====================" << std::endl;
  // Assigning Lambda expression for fourier sum to function gamma
  std::function<Vector2d(const double&)> gamma = [](const double& t){
    Vector2d res;
    res << cos(t) + 0.65*cos(2*t), 1.5*sin(t);
    return res;
  };
  // TODO: You may test your computation of the coefficients for gamma(t) and the
  // difference (gamma(0.1)-gamma(0))/||S_hat(0.1)-S_hat(0)||.
  Vector2d diffGamma; diffGamma.setZero();
  // Calculating gamma(s)-gamma(t) exactly
  Vector2d exDiffgamma = gamma(0.1) - gamma(0);
  std::cout << "(gamma (0.1) - gamma(0))/||S(0.1)-S(0)|| = " << std::endl;
  std::cout << "Calculated" << std::setw(20) << "Exact" <<  std::endl;
  std::cout << diffGamma.transpose() << std::setw(8)
      // Calculating (gamma(s)-gamma(t))/||S_hat(s)-S_hat(t)||
	    << exDiffgamma.transpose()/exDiffS.norm() << std::endl;
  std::cout << "============================================================="
	    << "=============" << std::endl << std::endl;



  //----------------------------------------------------------------------------
  std::cout << "================  Test L2-norm  ================"
	    << std::endl;
    VectorXd coeffToy(7);
  // Synthesized coefficients such that only cos(t) has a non zero coefficient
  coeffToy << 0,1,0,0,0,0,0;
  double l2norm = L2norm(Sprime, coeffToy);
  std::cout << "error computing L2 norm of cos(t) : "
	    << fabs( l2norm - std::sqrt(M_PI)) << std::endl
	    << "================================================"
	    << std::endl << std::endl;


  //----------------------------------------------------------------------------
  std::cout << "=====  Test system for gamma(t)  ====="
	    << std::endl;
  // Assigning Lambda expression to gammaprime according to gamma in 1.8.j
  std::function<Vector2d(const double&)> gammaprime = [](const double& t){
    Vector2d res;
    res << -sin(t) - 1.3*sin(2*t) , 1.5*cos(t);
    return res;
  };
  // Assigning Lambda expression to g corresponding to 1.8.j
  std::function<double(const Vector2d&)> g = [](const Vector2d& X){
    return sin(X(0))*sinh(X(1));
  };


  int Nl=15;
  VectorXi Nall(Nl);  Nall.setLinSpaced(Nl, 3, 31);
  VectorXd error(Nl); error.setZero();
  Vector2d T({0.5,0.2}); // Test point x in 1.8.j
    // TODO: Solve BIE for different Ns and compute error of the solution.

  // Output error and discretization parameters N
  std::ofstream out_error("AR_errors.txt");
  out_error << error;
  out_error.close( );

  std::ofstream out_N("AR_N.txt");
  out_N << Nall;
  out_N.close( );
  std::cout << "======================================"
	    << std::endl << std::endl;


  return 0;

}
