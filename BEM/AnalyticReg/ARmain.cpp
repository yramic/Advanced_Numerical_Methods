#include <iostream>
#include <fstream>
#include <istream> 
#include <cmath>
#include <Eigen/Dense>
// GMRES
#include <unsupported/Eigen/IterativeSolvers>

#if SOLUTION
#include "PeriodicTrapezoidalQR.hpp"
#endif // SOLUTION


/* SAM_LISTING_BEGIN_0 */
Eigen::MatrixXd computeAminusM(int N){
  Eigen::DiagonalMatrix<double, Eigen::Dynamic> AmM(2*N+1);

#if SOLUTION
  Eigen::VectorXd aux(2*N+1);
  aux << 0, Eigen::VectorXd::LinSpaced(N,1,N).cwiseInverse()*(M_PI/2.),
    Eigen::VectorXd::LinSpaced(N,1,N).cwiseInverse()*(M_PI/2.);
  AmM = aux.asDiagonal();

#else // TEMPLATE
  // TODO: Compute A-M
#endif // TEMPLATE
  return AmM;
}
/* SAM_LISTING_END_0 */


/* SAM_LISTING_BEGIN_1 */
  #if SOLUTION
template <typename PARAM>
Eigen::MatrixXd computeGammaCoefficients(const PARAM& gamma, int N){
  // Get quadrature points and weight
  Eigen::MatrixXd TR_points(2*N, 2);  double TR_w;
  std::tie(TR_points,TR_w) = PeriodicTrapRule(2*N);
  
  Eigen::MatrixXd coeffs(2*N,2);  coeffs.setZero();
  for(int k=1; k<=N; k++){
    // iterate over quadrature points
    for(int qp=0; qp<2*N; qp++){
      double z = TR_points(qp);
      Eigen::Vector2d gammaz = gamma(z);
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


Eigen::Vector2d evaluateGammaDiff(const Eigen::MatrixXd& gammaCoeffs, double s,
				  double t, int N){
  assert(gammaCoeffs.rows()==2*N);
  Eigen::VectorXd evalFun(2*N);  evalFun.setZero();
  // evaluate functions in series expansion
  for(int k=1; k<=N; k++){
    double Un; // this is actually U{n-1}
    if(fabs(s-t) < 1e-5)
      Un = k;
    else{
      Un = sin(k*(s-t)/2.)/(sin((s-t)/2.));
    }
    evalFun(k-1) = sin(k*(s+t)/2.)*Un;
    evalFun(k+N-1) = cos(k*(s+t)/2.)*Un;
  }

  // evaluate sum
  Eigen::Vector2d res = (-(gammaCoeffs.block(0,0,N,2)).transpose()
			 *evalFun.segment(0,N)
			 +(gammaCoeffs.block(N,0,N,2)).transpose()
			 *evalFun.segment(N,N) );
  return res;
}
#endif // SOLUTION

template <typename PARAM>
Eigen::MatrixXd computeM(const PARAM& gamma, int N){
  #if SOLUTION
  // Get quadrature points and weight
  int Nq = 8*N;
  Eigen::VectorXd TR_points(Nq);  double TR_w;
  std::tie(TR_points,TR_w) = PeriodicTrapRule(Nq);
  // For readibility, build each block separately and then assemble the
  // big matrix.
  Eigen::MatrixXd Mcc(N+1,N+1), Mss(N,N), Msc(N,N+1);
  Mcc.setZero(); Mss.setZero(); Msc.setZero(); 

  // In order to have a smooth integrand, we use the formula derived on 1.8.b
  // Get coefficients for Gamma in (1.0.42)
  Eigen::MatrixXd gammaCoeffs = computeGammaCoefficients(gamma,N);
  // Compute matrix entries
  for(int k=0; k<=N; k++){
    for(int l=0; l<=N; l++){
      // Compute double quadrature
      for(int qp1=0; qp1<Nq; qp1++){
	double s = TR_points(qp1);	  
	for(int qp2=0; qp2<Nq; qp2++){
	  double t = TR_points(qp2);
	  // compute argument of the log in the formula derived on 1.8.b
	  Eigen::Vector2d aux =  evaluateGammaDiff(gammaCoeffs, s, t, N);
	  // add contribution to the different matrices
	  Mcc(k,l)     += -TR_w/(4.*M_PI)*log(aux.squaredNorm())*cos(k*t)*cos(l*s);
	  
	  if(k>0){
	    Msc(k-1,l) += -TR_w/(4.*M_PI)*log(aux.squaredNorm())*sin(k*t)*cos(l*s);
	  }
	  
	  if(k>0 && l>0){
	    Mss(k-1,l-1) += -TR_w/(4.*M_PI)*log(aux.squaredNorm())*sin(k*t)*sin(l*s);
	  }
	  
	}// end loop over quadrature points
      }
    }// end for loop over l
  } // end for loop over k
  #endif // SOLUTION

  // Assemble big matrix
  Eigen::MatrixXd M(2*N+1, 2*N+1);

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
/* SAM_LISTING_END_1 */


/* SAM_LISTING_BEGIN_2 */
template <typename PARAM, typename FUNC>
Eigen::VectorXd computeG(const PARAM& gamma, const FUNC& g, int N){
  // Initialize right hand side vector
  Eigen::VectorXd RHS(2*N+1);  RHS.setZero();
  #if SOLUTION
  // Get quadrature points and weight
  Eigen::VectorXd TR_points(2*N);  double TR_w;
  std::tie(TR_points,TR_w) = PeriodicTrapRule(2*N);
  // Fill vector entries
  for(int k=0; k<=N; k++){
    // Evaluate integral by using the quadrature
    for(int qp=0; qp<2*N; qp++){
      double z = TR_points(qp);
      // First part (with $beta_k^c$)
      RHS(k) += TR_w * g(gamma(z))*cos(k*z);
      // Second part (with $beta_k^s$)
      if(k>0){
	RHS(k+N) += TR_w * g(gamma(z))*sin(k*z);
      }	
    }// end iteration over quadrature points
  }

  #else // TEMPLATE
    // TODO: Compute RHS
  #endif // TEMPLATE
  return RHS;
}
/* SAM_LISTING_END_2 */


/* SAM_LISTING_BEGIN_3 */
template <typename PARAM, typename FUNC>
Eigen::VectorXd solveBIE(const PARAM& gamma, const FUNC& g, int N){
    #if SOLUTION
  // In order to compute A we use A=(A-M)+M
  std::cout << " Assemble A " << std::endl; 
  Eigen::MatrixXd AmM = computeAminusM(N);
  Eigen::MatrixXd M   = computeM(gamma, N);
  // Build RHS
  std::cout << " Assemble RHS" << std::endl; 
  Eigen::VectorXd RHS = computeG(gamma, g, N);
  // Use direct solver
  std::cout << " Solve " << std::endl;
  Eigen::MatrixXd LHS = (AmM - M).eval();
  Eigen::VectorXd sol = LHS.lu().solve(RHS);
  std::cout << " Done " << std::endl;
#else // TEMPLATE
    // TODO: Build BIE system and solve it
  Eigen::VectorXd sol(2*N+1);
#endif // TEMPLATE
  
  return sol;
}
/* SAM_LISTING_END_3 */

  #if SOLUTION
template <typename PARAM, typename FUNC>
Eigen::VectorXd solveBIEonDisk(const PARAM& gamma, const FUNC& g, int N){
  // In order to compute A we do (A-M)+M
  std::cout << " Assemble A " << std::endl; 
  Eigen::MatrixXd AmM = computeAminusM(N);
  Eigen::MatrixXd M   = computeM(gamma, N);
  // Build RHS
  std::cout << " Assemble RHS" << std::endl; 
  Eigen::VectorXd RHS = computeG(gamma, g, N);
  std::cout << " Solve " << std::endl;
  Eigen::MatrixXd LHS = (AmM + M).block(1,1,2*N,2*N);
  Eigen::VectorXd sol = LHS.lu().solve(RHS.segment(1,2*N));
  std::cout << " Done " << std::endl;
  return sol;  
}



/* SAM_LISTING_BEGIN_4 */
template <typename PARAMDER>
double reconstructRho(const Eigen::VectorXd& coeffs, double t,
		      const PARAMDER& gammaprime){
  int N = (coeffs.rows()-1)/2; // asumming coeffs is a 2N+1 vector
  double res = coeffs(0);
  for(int k=1; k<=N; k++){
    res += coeffs(k)*cos(k*t)/(gammaprime(t)).norm()
      + coeffs(k+N)*sin(k*t)/(gammaprime(t)).norm();
  }
  return res;
}
#endif //SOLUTION
 
template <typename PARAMDER>
double L2norm(const PARAMDER& gammaprime, const Eigen::VectorXd& coeffs){
  double res=0.;
#if SOLUTION
  int N = (coeffs.rows()-1)/2.; // asumming coeffs is a 2N+1 vector
  // Get quadrature points and weight
  Eigen::VectorXd TR_points(2*N);  double TR_w;
  std::tie(TR_points,TR_w) = PeriodicTrapRule(4*N);
  for(int qp=0; qp<4*N; qp++){
    auto z = TR_points(qp);
    // evaluate
    double aux = reconstructRho(coeffs, z, gammaprime);
    res += TR_w*aux*aux;      
  }
  #else // TEMPLATE
    // TODO: reconstruct the function and compute its L2norm
  #endif // TEMPLATE
  
  return std::sqrt(res);
}
/* SAM_LISTING_END_4 */


/* SAM_LISTING_BEGIN_5 */
 #if SOLUTION
template <typename PARAM, typename PARAMDER>
double repFormulaSL(const Eigen::VectorXd& mu, const Eigen::Vector2d& X,
		    const PARAM& gamma,
		    const PARAMDER& gammaprime){
  int N = (mu.rows()-1)/2; // asumming coeffs is a 2N+1 vector
  double res=0.;
  Eigen::VectorXd TR_points(2*N);  double TR_w;
  std::tie(TR_points,TR_w) = PeriodicTrapRule(2*N);
  for(int qp=0; qp<2*N; qp++){
    auto z = TR_points(qp); 
    double aux = reconstructRho(mu, z, gammaprime);
    res += -TR_w/(2*M_PI)*log((X- gamma(z)).norm())*aux*(gammaprime(z)).norm();
  }
  return res;
}
 #endif //SOLUTION
/* SAM_LISTING_END_5 */

int main() {

  int N = 9;
  //----------------------------------------------------------------------------
  std::cout << "===========  Test integration  ===========" << std::endl;
  double Qint1 = 0., Qintcos = 0., Qintlogcos2 = 0.;
#if SOLUTION
  Eigen::MatrixXd TR_points(N, 2);  double TR_w;
  std::tie(TR_points,TR_w) = PeriodicTrapRule(N);
  for(int qp=0; qp<N; qp++){
    auto z = TR_points(qp);
    Qint1 += TR_w;
    Qintcos += TR_w*cos(z);
    Qintlogcos2 += TR_w*log(cos(z)*cos(z)+1.);
  }
#else // TEMPLATE
  // TODO: You may test your quadrature by integrating 1, cos(t) and log(cos(t)^2 + 1)
  #endif // TEMPLATE
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
  std::function<Eigen::Vector2d(const double&)> S = [](const double& t){
    Eigen::Vector2d res;
    res << cos(t) , sin(t);
    return res;
  };
#if SOLUTION
  Eigen::MatrixXd SCoeffs = computeGammaCoefficients(S,N);
  Eigen::Vector2d diffS = evaluateGammaDiff(SCoeffs, 0.1, 0, N);
#else // TEMPLATE
  // TODO: You may test your computation of the coefficients for S(t) and the
  // difference (S(0.1)-S(0))/S(0.1)-S(0).
  Eigen::Vector2d diffS; diffS.setZero();
  #endif // TEMPLATE
  Eigen::Vector2d exDiffS; exDiffS<< cos(0.1) - cos(0), sin(0.1);
  std::cout << "(S(0.1)-S(0))/S(0.1)-S(0) = " << diffS.transpose() << " vs "
	    << exDiffS(0)/exDiffS.norm() << " , " << exDiffS(1)/exDiffS.norm()
	    << std::endl;
  std::cout << "==============================================================" 
	    << "============" << std::endl << std::endl;


  //----------------------------------------------------------------------------
  std::cout << "====================  Test Coefficients for gamma(t)  "
	    << "====================" << std::endl;
  std::function<Eigen::Vector2d(const double&)> gamma = [](const double& t){
    Eigen::Vector2d res;
    res << cos(t) + 0.65*cos(2*t), 1.5*sin(t);
    return res;
  };
#if SOLUTION
  Eigen::MatrixXd gammaCoeffs = computeGammaCoefficients(gamma,N);
  Eigen::Vector2d diffGamma = evaluateGammaDiff(gammaCoeffs, 0.1, 0, N);
  #else // TEMPLATE
  // TODO: You may test your computation of the coefficients for S(t) and the
  // difference (S(0.1)-S(0))/S(0.1)-S(0).
  Eigen::Vector2d diffGamma; diffGamma.setZero();
  #endif // TEMPLATE
  Eigen::Vector2d exDiffgamma = gamma(0.1) - gamma(0);
  std::cout << "(gamma (0.1) - gamma(0))//S(0.1)-S(0) = " << diffGamma.transpose()
	    << " vs " << exDiffgamma.transpose()/exDiffS.norm() << std::endl;
  std::cout << "============================================================="
	    << "=============" << std::endl << std::endl;
  
#if SOLUTION
  //----------------------------------------------------------------------------
  std::cout << "=====  Test system for S(t) = (cos(t), sin(t))  ====="
	    << std::endl;
  std::function<Eigen::Vector2d(const double&)> Sprime = [](const double& t){
    Eigen::Vector2d res;
    res << -sin(t),  cos(t);
    return res;
  };

  std::function<double(const Eigen::Vector2d&)> gC = [](const Eigen::Vector2d& X){
    return X(0);
  };
  
  Eigen::VectorXd solC = solveBIEonDisk(S, gC, 10);
  Eigen::VectorXd solE(21); solE << 0, solC;
  double solCEval = reconstructRho(solE, M_PI, Sprime);
  std::cout << "eval at PI: " << solCEval
	    << " vs " << gC(S(M_PI))*2 << std::endl;
  std::cout << "====================================================="
	    << std::endl << std::endl;
  #endif // SOLUTION


  //----------------------------------------------------------------------------
  std::cout << "================  Test L2-norm  ================"
	    << std::endl;
    Eigen::VectorXd coeffToy(7);
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
  std::function<Eigen::Vector2d(const double&)> gammaprime = [](const double& t){
    Eigen::Vector2d res;
    res << -sin(t) - 1.3*sin(2*t) , 1.5*cos(t);
    return res;
  };
  std::function<double(const Eigen::Vector2d&)> g = [](const Eigen::Vector2d& X){
    return sin(X(0))*sinh(X(1));
  };

#if SOLUTION
  Eigen::VectorXi Nall(9); Nall<< 3,5,7,9,11,13,15,17,19;
  Eigen::VectorXd error(9); error.setZero();
  Eigen::Vector2d T({0.5,0.2});
  
  for(int j=0; j<9; j++){
    Eigen::VectorXd sol = solveBIE(gamma, g, Nall[j]);
    double solEval = repFormulaSL(sol, T, gamma, gammaprime);
    std::cout << j << ": " << solEval << " vs " << g(T)
	      << std::endl; 
    error(j) = fabs(solEval - g(T) );
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
