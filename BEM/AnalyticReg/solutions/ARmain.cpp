#include <iostream>
#include <fstream>
#include <istream> 
#include <cmath>
#include <Eigen/Dense>

using namespace Eigen;

#include "PeriodicTrapezoidalQR.hpp"


 
/* @brief Compute matrix A-M using analytic expression.
 * \param[in] N Discretization parameter indicating number of basis functions.
 */
/* SAM_LISTING_BEGIN_0 */
MatrixXd computeAminusM(int N){
  // Eigen offers a special type for diagonal matrices
  DiagonalMatrix<double, Dynamic> AmM(2*N+1);

  VectorXd aux(2*N+1);
  // Just the diagonal values from \eqref{arm:dm}
  aux << 0, VectorXd::LinSpaced(N,1,N).cwiseInverse()*(M_PI/2.),
    VectorXd::LinSpaced(N,1,N).cwiseInverse()*(M_PI/2.);
  // Convert vector into diagonal matrix
  AmM = aux.asDiagonal();

  return AmM;
}
/* SAM_LISTING_END_0 */


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
  for(int k=1; k<=N; k++){
    // iterate over quadrature points
    for(int qp=0; qp<2*N; qp++){
      double z = TR_points(qp);
      Vector2d gammaz = gamma(z);
      // compute coeficcients ak
      coeffs(k-1,0) += 1/M_PI*TR_w*gammaz(0)*cos(k*z);
      coeffs(k-1,1) += 1/M_PI*TR_w*gammaz(1)*cos(k*z);
      // compute coeficcients bk
      coeffs(k+N-1,0) += 1/M_PI*TR_w*gammaz(0)*sin(k*z);
      coeffs(k+N-1,1) += 1/M_PI*TR_w*gammaz(1)*sin(k*z);
    }// end for qp
  }// end for coeffs (k)
  return coeffs;  
}
/* SAM_LISTING_END_1a */


//----------------------------------------------------------------------------
/* @brief Evaluates the difference gamma(s)-gamma(t) according to its analytic 
 *        expression.
 * \param[in] gammaCoeffs Matrix containing 2N 2d-vectors corresponding to the 
 *                        Fourier expansion of gamma.
 * \param[in] s,t points in [0,2Pi]
 * \param[in] N Discretization parameter indicating number of basis functions.
 */
/* SAM_LISTING_BEGIN_1b */
Vector2d evaluateGammaDiff(const MatrixXd& gammaCoeffs,
			   double s, double t, int N){
  assert(gammaCoeffs.rows()==2*N);
  VectorXd evalFun(2*N); evalFun.setZero();
  // evaluate functions in series expansion
  for(int k=1; k<=N; k++){
    double Un; // this is actually $U_{n-1}$
    if(fabs(s-t) < 1e-5)  Un = k;
    else  Un = sin(k*(s-t)/2.)/(sin((s-t)/2.));

    evalFun(k-1) = sin(k*(s+t)/2.)*Un;
    evalFun(k+N-1) = cos(k*(s+t)/2.)*Un;
  }

  // evaluate sum
  Vector2d res = (-(gammaCoeffs.block(0,0,N,2)).transpose()
		  *evalFun.segment(0,N)
		  +(gammaCoeffs.block(N,0,N,2)).transpose()
		  *evalFun.segment(N,N));
  return res;
}
/* SAM_LISTING_END_1b */


//----------------------------------------------------------------------------
/* @brief Compute matrix M using periodic trapezoidal rule (2N+2 points).
 * \param[in] gamma Function that takes a double and returns a 2d vector  
 *                  corresponding to the parametrized curve.
 * \param[in] N Discretization parameter indicating number of basis functions.
 */
/* SAM_LISTING_BEGIN_1c */
template <typename PARAM>
MatrixXd computeM(const PARAM& gamma, int N){
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
      // Compute double quadrature
      for(int qp1=0; qp1<Nq; qp1++){
	double s = TR_points(qp1);	  
	for(int qp2=0; qp2<Nq; qp2++){
	  double t = TR_points(qp2);
	  // compute argument of the log in the formula derived on 1.8.b
	  Vector2d aux =  evaluateGammaDiff(gammaCoeffs, s, t, N);
	  // add contribution to the different matrices
	  Mcc(k,l)     += coeff*log(aux.squaredNorm())*cos(k*t)*cos(l*s);
	  
	  if(k>0){
	    Msc(k-1,l) += coeff*log(aux.squaredNorm())*sin(k*t)*cos(l*s);
	  }
	  
	  if(k>0 && l>0){
	    Mss(k-1,l-1) += coeff*log(aux.squaredNorm())*sin(k*t)*sin(l*s);
	  }
	  
	}// end loop over quadrature points
      }
    }// end for loop over l
  } // end for loop over k

  // Assemble big matrix
  MatrixXd M(2*N+1, 2*N+1);

  M.block(0  , 0  , N+1, N+1) = Mcc.eval();
  M.block(0  , N+1, N+1, N  ) = Msc.transpose().eval();
  M.block(N+1, 0  , N  , N+1) = Msc.eval();
  M.block(N+1, N+1, N  , N  ) = Mss.eval();

  return M;
}
/* SAM_LISTING_END_1c */


//----------------------------------------------------------------------------
/* @brief Compute right hand side using periodic trapezoidal rule (2N points).
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

  return RHS;
}
/* SAM_LISTING_END_2 */


//----------------------------------------------------------------------------
/* @brief Build and solve boundary integral equation V rho = g
 * \param[in] gamma Function that takes a double and returns a 2d vector  
 *                  corresponding to the parametrized curve.
 * \param[in] g Right hand side function (takes 2d points and returns a double).
 * \param[in] N Discretization parameter indicating number of basis functions.
 */
/* SAM_LISTING_BEGIN_3 */
template <typename PARAM, typename FUNC>
VectorXd solveBIE(const PARAM& gamma, const FUNC& g, int N){
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
  
  return sol;
}
/* SAM_LISTING_END_3 */


//----------------------------------------------------------------------------
/* @brief Build and solve boundary integral equation V rho = g ignoring first
 *        row of the system (A has empty row for on the Disk).
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
  int N = (coeffs.rows()-1)/2; // asumming coeffs is a 2N+1 vector
  double res = coeffs(0);
  for(int k=1; k<=N; k++){
    res += coeffs(k)*cos(k*t)/(gammaprime(t)).norm()
      + coeffs(k+N)*sin(k*t)/(gammaprime(t)).norm();
  }
  return res;
}
/* SAM_LISTING_END_4a */


//----------------------------------------------------------------------------
/* @brief Compute L2 norm of UN from its coefficients using periodic trapezoidal 
 *        rule (2N points).
 * \param[in] gammaprime Function that takes a double and returns a 2d vector  
 *                       corresponding to the derivative of the curve's 
 *                       parametrization.
 * \param[in] coeffs coefficients of UN
 */
/* SAM_LISTING_BEGIN_4b */
template <typename PARAMDER>
double L2norm(const PARAMDER& gammaprime, const VectorXd& coeffs){
  double res=0.;
  int N = (coeffs.rows()-1)/2.; // asumming coeffs is a 2N+1 vector
  // Get quadrature points and weight
  VectorXd TR_points(2*N);  double TR_w;
  std::tie(TR_points,TR_w) = PeriodicTrapRule(2*N);
  for(int qp=0; qp<2*N; qp++){
    auto z = TR_points(qp);
    // evaluate
    double aux = reconstructRho(coeffs, z, gammaprime);
    res += TR_w*aux*aux;      
  }
  
  return std::sqrt(res);
}
/* SAM_LISTING_END_4b */


//----------------------------------------------------------------------------
/* @brief Evaluate Single Layer Potential of function given by the coefficient 
 *                 vector mu on the point X and using the parametrization gamma.
 *                 The integration is done using periodic trapezoidal rule 
 *                 (2N+2 points).
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
  int N = (mu.rows()-1)/2; // asumming coeffs is a 2N+1 vector
  double res=0.;
  int Nq = 2*N+2;
  VectorXd TR_points(Nq);  double TR_w;
  std::tie(TR_points,TR_w) = PeriodicTrapRule(Nq);
  for(int qp=0; qp<Nq; qp++){
    auto z = TR_points(qp); 
    double aux = reconstructRho(mu, z, gammaprime);
    res += -TR_w/(2*M_PI)*log((X- gamma(z)).norm())*aux*(gammaprime(z)).norm();
  }
  return res;
}
/* SAM_LISTING_END_5 */



int main() {

  int N = 10;
  //----------------------------------------------------------------------------
  std::cout << "===========  Test integration  ===========" << std::endl;
  double Qint1 = 0., Qintcos = 0., Qintlogcos2 = 0.;
  MatrixXd TR_points(N, 2);  double TR_w;
  std::tie(TR_points,TR_w) = PeriodicTrapRule(N);
  for(int qp=0; qp<N; qp++){
    auto z = TR_points(qp);
    Qint1 += TR_w;
    Qintcos += TR_w*cos(z);
    Qintlogcos2 += TR_w*log(cos(z)*cos(z)+1.);
  }
  double intlogcos2ex = 2.*M_PI*log((sqrt(2)+1)/((sqrt(2)-1)*4));
  std::cout << " int 1                : " << fabs (Qint1 - 2*M_PI)
	    << "\n int cosx             : " << fabs(Qintcos)
	    << "\n int log(cosx*cosx+1) : " << fabs(Qintlogcos2 - intlogcos2ex)
	    << std::endl;
  std::cout << "==========================================" << std::endl
	    << std::endl;

  
  //----------------------------------------------------------------------------
  std::cout << "============  Test Coefficients for S(t) = (cos(t), sin(t))  "
	    << "============="  << std::endl;
  N = 4;
  std::function<Vector2d(const double&)> S = [](const double& t){
    Vector2d res;
    res << cos(t) , sin(t);
    return res;
  };
  std::function<Vector2d(const double&)> Sprime = [](const double& t){
    Vector2d res;
    res << -sin(t),  cos(t);
    return res;
  };

  
  MatrixXd SCoeffs = computeGammaCoefficients(S,N);
  Vector2d diffS = evaluateGammaDiff(SCoeffs, 0.1, 0, N);
  Vector2d exDiffS; exDiffS<< cos(0.1) - cos(0), sin(0.1);
  std::cout << "(S(0.1)-S(0))/S(0.1)-S(0) = " << diffS.transpose() << " vs "
	    << exDiffS(0)/exDiffS.norm() << " , " << exDiffS(1)/exDiffS.norm()
	    << std::endl;
  std::cout << "==============================================================" 
	    << "============" << std::endl << std::endl;


  //----------------------------------------------------------------------------
  std::cout << "====================  Test Coefficients for gamma(t)  "
	    << "====================" << std::endl;
  std::function<Vector2d(const double&)> gamma = [](const double& t){
    Vector2d res;
    res << cos(t) + 0.65*cos(2*t), 1.5*sin(t);
    return res;
  };
  MatrixXd gammaCoeffs = computeGammaCoefficients(gamma,N);
  Vector2d diffGamma = evaluateGammaDiff(gammaCoeffs, 0.1, 0, N);
  Vector2d exDiffgamma = gamma(0.1) - gamma(0);
  std::cout << "(gamma (0.1) - gamma(0))//S(0.1)-S(0) = " << diffGamma.transpose()
	    << " vs " << exDiffgamma.transpose()/exDiffS.norm() << std::endl;
  std::cout << "============================================================="
	    << "=============" << std::endl << std::endl;
  
  //----------------------------------------------------------------------------
  std::cout << "=====  Test system for S(t) = (cos(t), sin(t))  ====="
	    << std::endl;
  std::function<double(const Vector2d&)> gC = [](const Vector2d& X){
    return X(0);
  };
  
  VectorXd solC = solveBIEonDisk(S, gC, 10);
  VectorXd solE(21); solE << 0, solC;
  double solCEval = reconstructRho(solE, M_PI, Sprime);
  std::cout << "ERROR for evaluating at PI: " << fabs(solCEval- gC(S(M_PI))*2 )
	    << std::endl;
  std::cout << "====================================================="
	    << std::endl << std::endl;


  //----------------------------------------------------------------------------
  std::cout << "================  Test L2-norm  ================"
	    << std::endl;
    VectorXd coeffToy(7);
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
  std::function<Vector2d(const double&)> gammaprime = [](const double& t){
    Vector2d res;
    res << -sin(t) - 1.3*sin(2*t) , 1.5*cos(t);
    return res;
  };
  std::function<double(const Vector2d&)> g = [](const Vector2d& X){
    return sin(X(0))*sinh(X(1));
  };


  int Nl=15;
  VectorXi Nall(Nl);  Nall.setLinSpaced(Nl, 3, 31);
  VectorXd error(Nl); error.setZero();
  Vector2d T({0.5,0.2});
  for(int j=0; j<Nl; j++){
    VectorXd sol = solveBIE(gamma, g, Nall[j]);
    double solEval = repFormulaSL(sol, T, gamma, gammaprime);
    error(j) = fabs(solEval - g(T) );
    std::cout << "Error on level " << j << ": " << error(j) << std::endl; 
  }

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
