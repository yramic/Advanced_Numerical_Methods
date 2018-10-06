#include <iostream>
#include <fstream>
#include <istream>
#include <cmath>
#include <Eigen/Dense>
#include <iomanip>

using namespace Eigen;

#if SOLUTION
#include "PeriodicTrapezoidalQR.hpp"
#endif // SOLUTION



/* @brief Compute matrix A-M using analytic expression.
 * \param[in] N Discretization parameter indicating number of basis functions.
 */
/* SAM_LISTING_BEGIN_0 */
MatrixXd computeAminusM(int N){
  // A-M is a diagonal matrix of size 2N+1 X 2N+1
  DiagonalMatrix<double, Dynamic> AminusM(2*N+1);

#if SOLUTION
  // Diagonal for A-M constructed as a vector
  VectorXd AminusM_diagonal(2*N+1);
  // Just the diagonal values from \eqref{arm:dm}
  AminusM_diagonal << 0,
                      VectorXd::LinSpaced(N,1,N).cwiseInverse()*(M_PI/2.),
                      VectorXd::LinSpaced(N,1,N).cwiseInverse()*(M_PI/2.);
  // Convert vector into diagonal matrix
  AminusM = AminusM_diagonal.asDiagonal();

#else // TEMPLATE
  // TODO: Compute A-M
#endif // TEMPLATE
  return AminusM;
}
/* SAM_LISTING_END_0 */


  #if SOLUTION
//----------------------------------------------------------------------------
/* @brief Compute Fourier coefficients corresponding to gamma.
 * \param[in] gamma Function that takes a double and returns a 2d vector
 *                  corresponding to the parametrized curve.
 * \param[in] N Discretization parameter indicating number of basis functions.
 */
/* SAM_LISTING_BEGIN_1a */
template <typename PARAM>
MatrixXd computeGammaCoefficients(const PARAM& gamma, int N){
  // Get quadrature points and weights
  MatrixXd TR_points(2*N, 2);  double TR_w;
  std::tie(TR_points,TR_w) = PeriodicTrapRule(2*N);
  // Variable for returning result.
  MatrixXd coeffs(2*N,2);  coeffs.setZero();
  // Matrices for storing cosine and sine coefficients
  MatrixXd cos_coefficients(N,2); cos_coefficients.setZero();
  MatrixXd sin_coefficients(N,2); sin_coefficients.setZero();
  // Getting Fourier coefficients using integrals evaluated by quadrature
  for(int k=1; k<=N; k++){
    // iterate over quadrature points
    for(int qp=0; qp<2*N; qp++){
      double z = TR_points(qp);
      Vector2d gammaz = gamma(z);
      // compute coefficients ak for cosine terms
      cos_coefficients(k-1,0) += 1/M_PI*TR_w*gammaz(0)*cos(k*z);
      cos_coefficients(k-1,1) += 1/M_PI*TR_w*gammaz(1)*cos(k*z);
      // compute coefficients bk for sine terms
      sin_coefficients(k-1,0) += 1/M_PI*TR_w*gammaz(0)*sin(k*z);
      sin_coefficients(k-1,1) += 1/M_PI*TR_w*gammaz(1)*sin(k*z);
    }// end for qp
  }// end for coeffs (k)
  // Stacking the cosine and sine coefficients into the final matrix
  coeffs << cos_coefficients,sin_coefficients;
  return coeffs;
}
/* SAM_LISTING_END_1a */


//----------------------------------------------------------------------------
/* @brief Evaluates the difference
 *        \f$\frac{\gamma(s)-\gamma(t)}{||\hat{S}(s)-\hat{S}(t)||}\f$
 *        according to its analytic expression.
 * \param[in] gammaCoeffs Matrix containing 2N 2d-vectors corresponding to the
 *                        Fourier expansion of gamma.
 * \param[in] s,t points in [0,2Pi]
 * \param[in] N Discretization parameter indicating number of basis functions.
 */
/* SAM_LISTING_BEGIN_1b */
Vector2d evaluateGammaDiff(const MatrixXd& gammaCoeffs,
			                     double s, double t, int N) {
  assert(gammaCoeffs.rows()==2*N);
  // Vector for storing sin and cos values in first and second half respectively
  VectorXd sine_values(N); sine_values.setZero();
  VectorXd cosine_values(N); cosine_values.setZero();
  // evaluating sine and cosine functions in series expansion
  for(int k=1; k<=N; k++){
    double Un; // this is actually $U_{n-1}$
    if(fabs(s-t) < 1e-5)
      Un = k; // Using limit for s-t -> 0
    else
      Un = sin(k*(s-t)/2.)/(sin((s-t)/2.));

    sine_values(k-1) = sin(k*(s+t)/2.)*Un;
    cosine_values(k-1) = cos(k*(s+t)/2.)*Un;
  }

  MatrixXd sine_coefficients = gammaCoeffs.block(N,0,N,2);
  MatrixXd cosine_coefficients = gammaCoeffs.block(0,0,N,2);
  // evaluate sum using matrix multiplication ([2XN] X [NX1] = [2X1])
  Vector2d res =
      (-(cosine_coefficients).transpose()*sine_values
		  + sine_coefficients.transpose()*cosine_values);
  return res;
}
#endif // SOLUTION
/* SAM_LISTING_END_1b */


//----------------------------------------------------------------------------
/* @brief Compute matrix M using periodic trapezoidal rule (2N+2 points).
 * \tparam PARAM Type for gamma which supports evaluation operator: gamma(double)
 * \param[in] gamma Function that takes a double and returns a 2d vector
 *                  corresponding to the parametrized curve.
 * \param[in] N Discretization parameter indicating number of basis functions.
 */
/* SAM_LISTING_BEGIN_1c */
template <typename PARAM>
MatrixXd computeM(const PARAM& gamma, int N){
  #if SOLUTION
  // Get quadrature points and weight
  int Nq = 2*N+2;
  VectorXd TR_points(Nq);  double TR_w;
  std::tie(TR_points,TR_w) = PeriodicTrapRule(Nq);
  // For readibility, build each block separately  and then assemble
  // the big matrix.
  MatrixXd Mcc(N+1,N+1), Mss(N,N), Msc(N,N+1);
  Mcc.setZero(); Mss.setZero(); Msc.setZero();

  // In order to have a smooth integrand, we use the formula derived on 1.8.b
  // Get coefficients for Gamma in (1.0.42)
  MatrixXd gammaCoeffs = computeGammaCoefficients(gamma,N);
  // Create variable multiplying quadrature weights and integral scaling
  double coeff = -TR_w*TR_w/(4.*M_PI);
  // Compute matrix entries
  for(int k=0; k<=N; k++){
    for(int l=0; l<=N; l++){
      // Computing entry using double quadrature
      for(int qp1=0; qp1<Nq; qp1++){
      	double s = TR_points(qp1);
      	for(int qp2=0; qp2<Nq; qp2++){
      	  double t = TR_points(qp2);
      	  // compute argument of the log in the formula derived on 1.8.b
      	  double arg =  evaluateGammaDiff(gammaCoeffs, s, t, N).squaredNorm();
      	  // add contribution to the different matrices
      	  Mcc(k,l)     += coeff*log(arg)*cos(k*t)*cos(l*s);

      	  if(k>0)
      	    Msc(k-1,l) += coeff*log(arg)*sin(k*t)*cos(l*s);

      	  if(k>0 && l>0)
      	    Mss(k-1,l-1) += coeff*log(arg)*sin(k*t)*sin(l*s);

      	}// end loop over quadrature points
      }
    }// end for loop over l
  } // end for loop over k
  #endif // SOLUTION

  // Assemble big matrix
  MatrixXd M(2*N+1, 2*N+1);

  #if SOLUTION
  M.block(0  , 0  , N+1, N+1) = Mcc.eval();
  M.block(0  , N+1, N+1, N  ) = Msc.transpose().eval();
  M.block(N+1, 0  , N  , N+1) = Msc.eval();
  M.block(N+1, N+1, N  , N  ) = Mss.eval();
  #else // TEMPLATE
    // TODO: Compute M using quadrature
  #endif // TEMPLATE

  return M;
}
/* SAM_LISTING_END_1c */


//----------------------------------------------------------------------------
/* @brief Compute right hand side using periodic trapezoidal rule (2N points).
 * \tparam PARAM Type for gamma which supports evaluation operator: gamma(double)
 * \tparam FUNC Type for function g which supports evaluation operator: g(double)
 * \param[in] gamma Function that takes a double and returns a 2d vector
 *                  corresponding to the parametrized curve.
 * \param[in] g Right hand side function (takes 2d points and returns a double).
 * \param[in] N Discretization parameter indicating number of basis functions.
 */
/* SAM_LISTING_BEGIN_2 */
template <typename PARAM, typename FUNC>
VectorXd computeG(const PARAM& gamma, const FUNC& g, int N){
  // Initialize right hand side vector
  VectorXd RHS(2*N+1);  RHS.setZero();
  #if SOLUTION
  // Get quadrature points and weight
  VectorXd TR_points(2*N);  double TR_w;
  std::tie(TR_points,TR_w) = PeriodicTrapRule(2*N);
  // Fill vector entries
  for(int k=0; k<=N; k++){
    // Evaluate integral by using the quadrature
    for(int qp=0; qp<2*N; qp++){
      double z = TR_points(qp);
      // First part (with $beta_k^c$)
      RHS(k) += TR_w * g(gamma(z))*cos(k*z);
      // Second part (with $beta_k^s$)
      if(k>0) RHS(k+N) += TR_w * g(gamma(z))*sin(k*z);
    }// end iteration over quadrature points
  }

  #else // TEMPLATE
    // TODO: Compute RHS
  #endif // TEMPLATE
  return RHS;
}
/* SAM_LISTING_END_2 */


//----------------------------------------------------------------------------
/* @brief Build and solve boundary integral equation V rho = g
 * \tparam PARAM Type for gamma which supports evaluation operator: gamma(double)
 * \tparam FUNC Type for function g which supports evaluation operator: g(double)
 * \param[in] gamma Function that takes a double and returns a 2d vector
 *                  corresponding to the parametrized curve.
 * \param[in] g Right hand side function (takes 2d points and returns a double).
 * \param[in] N Discretization parameter indicating number of basis functions.
 */
/* SAM_LISTING_BEGIN_3 */
template <typename PARAM, typename FUNC>
VectorXd solveBIE(const PARAM& gamma, const FUNC& g, int N){
    #if SOLUTION
  // In order to compute A we use A=(A-M)+M
  std::cout << " Assemble A " << std::endl;
  MatrixXd AmM = computeAminusM(N);
  MatrixXd M   = computeM(gamma, N);
  // Build RHS
  std::cout << " Assemble RHS" << std::endl;
  VectorXd RHS = computeG(gamma, g, N);
  // Use direct solver
  std::cout << " Solve " << std::endl;
  MatrixXd LHS = (AmM + M).eval();
  VectorXd sol = LHS.lu().solve(RHS);
  std::cout << " Done " << std::endl;
#else // TEMPLATE
    // TODO: Build BIE system and solve it
  VectorXd sol(2*N+1);
#endif // TEMPLATE

  return sol;
}
/* SAM_LISTING_END_3 */


//----------------------------------------------------------------------------
  #if SOLUTION
/* @brief Build and solve boundary integral equation V rho = g ignoring first
 *        row of the system (A has empty row for on the Disk).
 * \tparam PARAM Type for gamma which supports evaluation operator: gamma(double)
 * \tparam FUNC Type for function g which supports evaluation operator: g(double)
 * \param[in] gamma Function that takes a double and returns a 2d vector
 *                  corresponding to the parametrized curve.
 * \param[in] g Right hand side function (takes 2d points and returns a double).
 * \param[in] N Discretization parameter indicating number of basis functions.
 */
template <typename PARAM, typename FUNC>
VectorXd solveBIEonDisk(const PARAM& gamma, const FUNC& g, int N){
  // In order to compute A we do (A-M)+M
  std::cout << " Assemble A " << std::endl;
  MatrixXd AmM = computeAminusM(N);
  MatrixXd M   = computeM(gamma, N);
  // Build RHS
  std::cout << " Assemble RHS" << std::endl;
  VectorXd RHS = computeG(gamma, g, N);
  std::cout << " Solve " << std::endl;
  MatrixXd LHS = (AmM + M).block(1,1,2*N,2*N);
  VectorXd sol = LHS.lu().solve(RHS.segment(1,2*N));
  std::cout << " Done " << std::endl;
  return sol;
}


//----------------------------------------------------------------------------
/* @brief Reconstruct function UN from its coefficients and evaluate it at t.
 * \tparam PARAMDER Type for gammaprime which supports evaluation operator:
 *                  gammaprime(double)
 * \param[in] coeffs coefficients of UN
 * \param[in] t point in [0,2Pi]
 * \param[in] gammaprime Function that takes a double and returns a 2d vector
 *                       corresponding to the derivative of the curve's
 *                       parametrization.
 */
/* SAM_LISTING_BEGIN_4a */
template <typename PARAMDER>
double reconstructRho(const VectorXd& coeffs, double t,
		      const PARAMDER& gammaprime){
  int N = (coeffs.rows()-1)/2; // assumming coeffs is a 2N+1 vector
  double res = coeffs(0);
  for(int k=1; k<=N; k++){
    res += coeffs(k)*cos(k*t)/(gammaprime(t)).norm()
      + coeffs(k+N)*sin(k*t)/(gammaprime(t)).norm();
  }
  return res;
}
#endif //SOLUTION
/* SAM_LISTING_END_4a */


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
/* SAM_LISTING_BEGIN_4b */
template <typename PARAMDER>
double L2norm(const PARAMDER& gammaprime, const VectorXd& coeffs){
  double res=0.;
#if SOLUTION
  int N = (coeffs.rows()-1)/2.; // assumming coeffs is a 2N+1 vector
  // Get quadrature points and nodes
  VectorXd TR_points(2*N);  double TR_w;
  std::tie(TR_points,TR_w) = PeriodicTrapRule(2*N);
  for(int qp=0; qp<2*N; qp++){
    auto z = TR_points(qp);
    // evaluate
    double rho = reconstructRho(coeffs, z, gammaprime);
    res += TR_w*rho*rho;
  }
  #else // TEMPLATE
    // TODO: reconstruct the function and compute its L2norm
  #endif // TEMPLATE

  return std::sqrt(res);
}
/* SAM_LISTING_END_4b */


 #if SOLUTION
//----------------------------------------------------------------------------
/* @brief Evaluate Single Layer Potential of function given by the coefficient
 *                 vector mu on the point X and using the parametrization gamma.
 *                 The integration is done using periodic trapezoidal rule
 *                 (2N+2 points).
 * \tparam PARAM Type for gamma which supports evaluation operator: gamma(double)
 * \tparam PARAMDER Type for gammaprime which supports evaluation operator:
 *                  gammaprime(double)
 * \param[in] coeffs coefficients of UN.
 * \param[in] X 2d vector containing evaluation point.
 * \param[in] gamma Function that takes a double and returns a 2d vector
 *                  corresponding to the parametrized curve.
 * \param[in] gammaprime Function that takes a double and returns a 2d vector
 *                       corresponding to the derivative of the curve's
 *                       parametrization.
 */
/* SAM_LISTING_BEGIN_5 */
template <typename PARAM, typename PARAMDER>
double repFormulaSL(const VectorXd& mu, const Vector2d& X,
		    const PARAM& gamma,
		    const PARAMDER& gammaprime){
  int N = (mu.rows()-1)/2; // assumming coeffs is a 2N+1 vector
  double res=0.; // For storing the constructed single layer potential
  int Nq = 2*N+2;
  // Getting quadrature weights and nodes for periodic trapezoidal rule
  VectorXd TR_points(Nq);  double TR_w;
  std::tie(TR_points,TR_w) = PeriodicTrapRule(Nq);
  // Performing integration for single layer potential via quadrature
  for(int qp=0; qp<Nq; qp++){
    auto z = TR_points(qp);
    double rho = reconstructRho(mu, z, gammaprime);
    // The single layer potential formula
    res += -TR_w/(2*M_PI)*log((X- gamma(z)).norm())*rho*(gammaprime(z)).norm();
  }
  return res;
}
 #endif //SOLUTION
/* SAM_LISTING_END_5 */



int main() {

  int N = 10;
  //----------------------------------------------------------------------------
  std::cout << "===========  Test integration errors ===========" << std::endl;
  double Qint1 = 0., Qintcos = 0., Qintlogcos2 = 0.;
#if SOLUTION
  // Getting quadrature weights and nodes for periodic trapezoidal rule
  MatrixXd TR_points(N, 2);  double TR_w;
  std::tie(TR_points,TR_w) = PeriodicTrapRule(N);
  // Performing integration via quadrature
  for(int qp=0; qp<N; qp++){
    auto z = TR_points(qp);
    Qint1 += TR_w;
    Qintcos += TR_w*cos(z);
    Qintlogcos2 += TR_w*log(cos(z)*cos(z)+1.);
  }
#else // TEMPLATE
  /* TODO: You may test your quadrature by integrating 1, cos(t) and
     log(cos(t)^2 + 1) */
  #endif // TEMPLATE
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


#if SOLUTION
  MatrixXd SCoeffs = computeGammaCoefficients(S,N); // gamma = S = S_hat
  Vector2d diffS = evaluateGammaDiff(SCoeffs, 0.1, 0, N);
#else // TEMPLATE
  // TODO: You may test your computation of the coefficients for gamma(t) and the
  // difference (gamma(0.1)-gamma(0))/||S_hat(0.1)-S_hat(0)||.
  Vector2d diffS; diffS.setZero();
#endif // TEMPLATE
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
#if SOLUTION
  MatrixXd gammaCoeffs = computeGammaCoefficients(gamma,N);
  Vector2d diffGamma = evaluateGammaDiff(gammaCoeffs, 0.1, 0, N);
  #else // TEMPLATE
  // TODO: You may test your computation of the coefficients for gamma(t) and the
  // difference (gamma(0.1)-gamma(0))/||S_hat(0.1)-S_hat(0)||.
  Vector2d diffGamma; diffGamma.setZero();
  #endif // TEMPLATE
  // Calculating gamma(s)-gamma(t) exactly
  Vector2d exDiffgamma = gamma(0.1) - gamma(0);
  std::cout << "(gamma (0.1) - gamma(0))/||S(0.1)-S(0)|| = " << std::endl;
  std::cout << "Calculated" << std::setw(20) << "Exact" <<  std::endl;
  std::cout << diffGamma.transpose() << std::setw(8)
      // Calculating (gamma(s)-gamma(t))/||S_hat(s)-S_hat(t)||
	    << exDiffgamma.transpose()/exDiffS.norm() << std::endl;
  std::cout << "============================================================="
	    << "=============" << std::endl << std::endl;

#if SOLUTION
  //----------------------------------------------------------------------------
  std::cout << "=====  Test system for S(t) = (cos(t), sin(t))  ====="
	    << std::endl;
  // Assigning Lambda expression to the function gC
  std::function<double(const Vector2d&)> gC = [](const Vector2d& X){
    return X(0);
  };

  // Getting coefficients
  VectorXd solC = solveBIEonDisk(S, gC, 10);
  // Adding 0 to the solution as solveBIEonDisk ignores the empty first row
  VectorXd solE(21); solE << 0, solC;
  // Reconstructing rho from the coefficients
  double solCEval = reconstructRho(solE, M_PI, Sprime);
  std::cout << "ERROR for evaluating at PI: " << fabs(solCEval- gC(S(M_PI))*2 )
	    << std::endl;
  std::cout << "====================================================="
	    << std::endl << std::endl;
  #endif // SOLUTION


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
  /* SAM_LISTING_BEGIN_6 */
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
#if SOLUTION
  for(int j=0; j<Nl; j++){
    // Solving BIE to get coefficients
    VectorXd sol = solveBIE(gamma, g, Nall[j]);
    // Evaluating Single Layer potential at point T using obtained coefficients
    double solEval = repFormulaSL(sol, T, gamma, gammaprime);
    error(j) = fabs(solEval - g(T) );
    std::cout << "Error on level " << j << ": " << error(j) << std::endl;
  }
#else // TEMPLATE
    // TODO: Solve BIE for different Ns and compute error of the solution.
#endif // TEMPLATE

  // Output error and discretization parameters N
  std::ofstream out_error("AR_errors.txt");
  out_error << error;
  out_error.close( );

  std::ofstream out_N("AR_N.txt");
  out_N << Nall;
  out_N.close( );
  std::cout << "======================================"
	    << std::endl << std::endl;

  /* SAM_LISTING_END_6 */

  return 0;

}
