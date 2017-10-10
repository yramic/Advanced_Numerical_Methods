#include <iostream>
#include <cmath>
#include <Eigen/Dense>
// GMRES
#include <unsupported/Eigen/IterativeSolvers>
#include "PeriodicTrapezoidalQR.hpp"


/* SAM_LISTING_BEGIN_0 */
Eigen::MatrixXd computeAminusM(int N){
  Eigen::DiagonalMatrix<double, Eigen::Dynamic> AmM(2*N+1);
  Eigen::VectorXd aux(2*N+1);
  aux << 0, Eigen::VectorXd::LinSpaced(N,1,N).cwiseInverse()*(-M_PI/2.),
    Eigen::VectorXd::LinSpaced(N,1,N).cwiseInverse()*(-M_PI/2.);
  AmM = aux.asDiagonal();
  return AmM;
}
/* SAM_LISTING_END_0 */


/* SAM_LISTING_BEGIN_1 */
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

Eigen::VectorXd evaluateCosSinUn(double t, double s, int N){
  Eigen::VectorXd res(2*N);  res.setZero();
  for(int k=1; k<=N; k++){
    double Un; // this is actually U_{n-1}
    if(fabs(s-t) < 1e-5)
      Un = k;
    else{
      Un = sin(k*(t-s)/2.)/(sin((t-s)/2.));
    }
    res(k-1) = sin(k*(s+t)/2.)*Un;
    res(k+N-1) = cos(k*(s+t)/2.)*Un;
  }
  return res;
}

template <typename PARAM>
Eigen::MatrixXd computeM(const PARAM& gamma, int N){
  // Get quadrature points and weight
  int Nq = 4*N;
  Eigen::VectorXd TR_points(Nq);  double TR_w;
  std::tie(TR_points,TR_w) = PeriodicTrapRule(Nq);
  // For readibility, build each separate block and then assemble the big matrix.
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
	  // compute argument of the log in hat{k}
	  Eigen::VectorXd evalCS = evaluateCosSinUn(t,s, N);
	  Eigen::Vector2d aux = (-(gammaCoeffs.block(0,0,N,2)).transpose()
				 *evalCS.segment(0,N)
		                 +(gammaCoeffs.block(N,0,N,2)).transpose()
				 *evalCS.segment(N,N) );

	  // add contribution to the different matrices
	  double valcc = -TR_w/(4.*M_PI)*log(aux.squaredNorm())*cos(k*t)*cos(l*s);
	  double valsc = 0.; double valss =0.;
	  Mcc(k,l)     += valcc;
	  
	  if(k>0){
	    valsc = -TR_w/(4.*M_PI)*log(aux.squaredNorm())*sin(k*t)*cos(l*s);
	    Msc(k-1,l) += valsc;
	  }
	  
	  if(k>0 && l>0){
	    valss = -TR_w/(4.*M_PI)*log(aux.squaredNorm())*sin(k*t)*sin(l*s);
	    Mss(k-1,l-1) += valss;
	  }
	  
	}// end loop over quadrature points
      }
    }// end for loop over l
  } // end for loop over k

  // Assemble big matrix
  Eigen::MatrixXd M(2*N+1, 2*N+1);
  M.block(0  , 0  , N+1, N+1) = Mcc.eval();
  M.block(0  , N+1, N+1, N  ) = Msc.transpose().eval();
  M.block(N+1, 0  , N  , N+1) = Msc.eval();
  M.block(N+1, N+1, N  , N  ) = Mss.eval();

  return M;
}
/* SAM_LISTING_END_1 */


/* SAM_LISTING_BEGIN_2 */
template <typename PARAM, typename FUNC>
Eigen::VectorXd computeG(const PARAM& gamma, const FUNC& g, int N){
  // Initialize right hand side vector
  Eigen::VectorXd RHS(2*N+1);  RHS.setZero();
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
  return RHS;
}
/* SAM_LISTING_END_2 */


/* SAM_LISTING_BEGIN_3 */
template <typename PARAM, typename FUNC>
Eigen::VectorXd solveBIE(const PARAM& gamma, const FUNC& g, int N){
  // In order to compute A we do (A-M)+M
  std::cout << " Assemble A " << std::endl; 
  Eigen::MatrixXd AmM = computeAminusM(N);
  Eigen::MatrixXd M   = computeM(gamma, N);
  std::cout << "M norm : " << M.norm() << " , " << M.cols() <<std::endl;
  std::cout << "A-M norm : " << AmM.norm() << " , " << M.cols() <<std::endl;
  // Build RHS
  std::cout << " Assemble RHS" << std::endl; 
  Eigen::VectorXd RHS = computeG(gamma, g, N);
  std::cout << " Solve " << std::endl;
  // Declare GMRES solver
  Eigen::ConjugateGradient< Eigen::MatrixXd > solver;
  solver.setMaxIterations( M.cols() );
  solver.setTolerance( 1.e-5 );
  // initialize the solver
  Eigen::MatrixXd LHS = (AmM + M).eval();
  solver.compute( LHS );
  Eigen::VectorXd sol = solver.solve(RHS);
    std::cout << "  #iterations   = " << solver.iterations() << std::endl
	      << "  estimated error = " << solver.error()    << std::endl;
  std::cout << " Done " << std::endl; 
  return sol;
}

template <typename PARAM, typename FUNC>
Eigen::VectorXd solveBIE2(const PARAM& gamma, const FUNC& g, int N){
  // In order to compute A we do (A-M)+M
  std::cout << " Assemble A " << std::endl; 
  Eigen::MatrixXd AmM = computeAminusM(N);
  Eigen::MatrixXd M   = computeM(gamma, N);
  std::cout << "M norm : " << M.norm() << " , " << M.cols() <<std::endl;
  std::cout << "A-M norm : " << AmM.norm() << " , " << M.cols() <<std::endl;
  // Build RHS
  std::cout << " Assemble RHS" << std::endl; 
  Eigen::VectorXd RHS = computeG(gamma, g, N);
  std::cout << " Solve " << std::endl;
  // Declare GMRES solver
  Eigen::ConjugateGradient< Eigen::MatrixXd > solver;
  solver.setMaxIterations( M.cols() );
  solver.setTolerance( 1.e-5 );
  // initialize the solver
  Eigen::MatrixXd LHS = (AmM + M).block(1,1,2*N,2*N);
  solver.compute( LHS );
  Eigen::VectorXd sol = solver.solve(RHS.segment(1,2*N));
    std::cout << "  #iterations   = " << solver.iterations() << std::endl
	      << "  estimated error = " << solver.error()    << std::endl;
  std::cout << " Done " << std::endl; 
  return sol;
}
/* SAM_LISTING_END_3 */


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

template <typename PARAMDER>
double L2norm(const PARAMDER& gammaprime, const Eigen::VectorXd& coeffs){
  double res=0.;
  int N = (coeffs.rows()-1)/2.; // asumming coeffs is a 2N+1 vector
  // Get quadrature points and weight
  Eigen::VectorXd TR_points(4*N);  double TR_w;
  std::tie(TR_points,TR_w) = PeriodicTrapRule(4*N);
  for(int qp=0; qp<4*N; qp++){
    auto z = TR_points(qp);
    // evaluate
    double aux = reconstructRho(coeffs, z, gammaprime);
    res += TR_w*aux*aux;      
  }
  return std::sqrt(res);
}
/* SAM_LISTING_END_4 */

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

int main() {

 std::cout << "==========  Test integration  ==========" << std::endl;
 int N = 5;
  Eigen::MatrixXd TR_points(N, 2);  double TR_w;
  std::tie(TR_points,TR_w) = PeriodicTrapRule(N);
  double Qint1 = 0.;  double Qintx = 0.;  double Qintcos = 0.;
  double Qintlogcos2 = 0.;
  for(int qp=0; qp<N; qp++){
    auto z = TR_points(qp);
    Qint1 += TR_w;
    Qintcos += TR_w*cos(z);
    Qintlogcos2 += TR_w*log(cos(z)*cos(z)+2.);
  }
  std::cout << " int 1              : " << Qint1 << " vs " << 2*M_PI
	    << "\n int cosx           :"  << Qintcos
	    << " vs 0 \n int log(cosx*cosx) : " << Qintlogcos2 << " vs "
	    << 5.693428 << std::endl;
  std::cout << "========================================" << std::endl << std::endl;

  
  std::cout << "============  Test Coefficients for S(t) = (cos(t), sin(t))  ============="
	    << std::endl;
    std::function<Eigen::Vector2d(const double&)> S = [](const double& t){
    Eigen::Vector2d res;
    res << cos(t) , sin(t);
    return res;
  };
    
  Eigen::MatrixXd SCoeffs = computeGammaCoefficients(S,N);
  Eigen::VectorXd evalFun = evaluateCosSinUn(0.1, 0, N);
  Eigen::Vector2d diffS = (-SCoeffs.block(0,0,N,2)).transpose()*evalFun.segment(0,N)
              +(SCoeffs.block(N,0,N,2)).transpose()*evalFun.segment(N,N);
  Eigen::Vector2d exDiffS; exDiffS<< cos(0.1) - cos(0), sin(0.1);
  std::cout << "(S(0.1)-S(0))/S(0.1)-S(0) = " << diffS.transpose() << " vs "
	    << exDiffS(0)/exDiffS.norm() << " , " << exDiffS(1)/exDiffS.norm()
	    << std::endl;
  std::cout << "=========================================================================="
	    << std::endl << std::endl;

  
  std::cout << "=============  Test Coefficients for gamma as in assignment  ============="
	    << std::endl;
  std::function<Eigen::Vector2d(const double&)> gamma = [](const double& t){
    Eigen::Vector2d res;
    res << cos(t) + 0.65*cos(2*t), 1.5*sin(t);
    return res;
  };

  Eigen::MatrixXd gammaCoeffs = computeGammaCoefficients(gamma,N);
  Eigen::Vector2d diffGamma = (-gammaCoeffs.block(0,0,N,2)).transpose()
                               *evalFun.segment(0,N)
                      +(gammaCoeffs.block(N,0,N,2)).transpose()*evalFun.segment(N,N);
  Eigen::Vector2d exDiffgamma = gamma(0.1) - gamma(0);
  std::cout << "(gamma (0.1) - gamma(0))//S(0.1)-S(0) = " << diffGamma.transpose()
	    << " vs " << exDiffgamma.transpose()/exDiffS.norm() << std::endl;
  std::cout << "=========================================================================="
	    << std::endl << std::endl;
  

  std::cout << "=============  Test system for S(t) = (cos(t), sin(t))  ============="
	     << std::endl;
  std::function<Eigen::Vector2d(const double&)> Sprime = [](const double& t){
    Eigen::Vector2d res;
    res << -sin(t),  cos(t);
    return res;
  };

  std::function<double(const Eigen::Vector2d&)> gC = [](const Eigen::Vector2d& X){
    return X(0);
  };
  
  Eigen::VectorXd solC = solveBIE2(S, gC, 10);
  double solCEval = reconstructRho(solC, M_PI, Sprime);
  std::cout << "eval at PI: " << solCEval
	    << " vs " << gC(S(M_PI))*2 << std::endl;
  std::cout << "====================================================================="
	    << std::endl << std::endl;

  std::cout << "=============  Test system for S2(t) = (2cos(t), 2sin(t))  ============="
	    << std::endl;
   std::function<Eigen::Vector2d(const double&)> S2 = [](const double& t){
    Eigen::Vector2d res;
    res << 2*cos(t), 2*sin(t);
    return res;
  };
  std::function<Eigen::Vector2d(const double&)> S2prime = [](const double& t){
    Eigen::Vector2d res;
    res << -2*sin(t),  2*cos(t);
    return res;
  };
  
  Eigen::VectorXd solC2 = solveBIE2(S2, gC, 10);
  double solC2Eval = reconstructRho(solC2, M_PI, S2prime);
  std::cout << "eval at PI: " << solC2Eval
	    << " vs " << gC(S2(M_PI))*2 << std::endl;
  std::cout << "====================================================================="
	    << std::endl << std::endl;
  
  
  std::cout << "================  Test L2-norm  ================"
	    << std::endl;
    Eigen::VectorXd coeffToy(7);
  coeffToy << 0,1,0,0,0,0,0;
  double l2norm = L2norm(Sprime, coeffToy);
  std::cout << "L2 norm of cos(t) : " << l2norm << " vs " << std::sqrt(M_PI)
	    << std::endl
	    << "================================================"
	    << std::endl << std::endl;
    
  
  // SYSTEM
  std::function<Eigen::Vector2d(const double&)> gammaprime = [](const double& t){
    Eigen::Vector2d res;
    res << -sin(t) - 1.3*sin(2*t) , 1.5*cos(t);
    return res;
  };

  std::function<double(const Eigen::Vector2d&)> g = [](const Eigen::Vector2d& X){
    return sin(X(0))*sinh(X(1));
  };


  std::vector<int> Nall({3,5,7,9,11});
  int NumN = Nall.size();
  Eigen::VectorXd error(NumN);
  Eigen::Vector2d T({0.5,0.2});
  for(int j=0; j<NumN; j++){
    Eigen::VectorXd sol = solveBIE(gamma, g, Nall[j]);
    double solEval = repFormulaSL(sol, T, gamma, gammaprime);
    std::cout << j << ": " << solEval << " vs " << g(T)
	      << std::endl; 
    error(j) = fabs(solEval + g(T) );
  }
  
  std::cout << error << std::endl;  
  
    
  return 0;

}
