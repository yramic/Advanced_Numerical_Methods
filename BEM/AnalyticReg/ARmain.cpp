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
  aux << 0, Eigen::VectorXd::LinSpaced(N,1,N).cwiseInverse()*(M_PI/2.),
    Eigen::VectorXd::LinSpaced(N,1,N).cwiseInverse()*(M_PI/2.);
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
  
  Eigen::MatrixXd coeffs(2*N+1,2);  coeffs.setZero();
  for(int k=0; k<=N; k++){
    // iterate over quadrature points
    for(int qp=0; qp<2*N; qp++){
      auto z = TR_points(qp);
      auto gammaz = gamma(z);
      // compute coeficcients ak
      coeffs(k,0) += 1/M_PI*TR_w*gammaz(0)*cos(k*z);
      coeffs(k,1) += 1/M_PI*TR_w*gammaz(1)*cos(k*z);

      // compute coefficients bk
      if(k>0){
	coeffs(k+N,0) += 1/M_PI*TR_w*gammaz(0)*sin(k*z);
	coeffs(k+N,1) += 1/M_PI*TR_w*gammaz(1)*sin(k*z);
      }
    }// end for qp
  }// end for coeffs (k)
  return coeffs;  
}

Eigen::VectorXd evaluateCosSinUn(double t, double s, int N){
  Eigen::VectorXd res(2*N+1);  res.setZero();
  for(int k=0; k<=N; k++){
    double Un;
    if(fabs(s-t) < 1e-3)
      Un = 1.+k;
    else{
      if(k==0) 
	Un = 0;
      else
	Un = sin(k*(t-s)/2.)/(sin((t-s)/2.));
    }

    res(k) = sin(k*(s+t)/2.)*Un;
    if(k>0)
      res(k+N) = cos(k*(s+t)/2.)*Un;
  }
  return res;
}


template <typename PARAM>
Eigen::MatrixXd computeM(const PARAM& gamma, int N){
  // Get quadrature points and weight
  Eigen::VectorXd TR_points(4*N);  double TR_w;
  std::tie(TR_points,TR_w) = PeriodicTrapRule(4*N);
  // For readibility, build each separate block and then assemble the big matrix.
  Eigen::MatrixXd Mcc(N+1,N+1), Mss(N,N), Msc(N,N+1);
  Mcc.setZero(); Mss.setZero(); Msc.setZero(); 

  // In order to have a smooth integrand, we use the formula derived on 1.8.b
  // Get coefficients for Gamma in (1.0.42)
  auto gammaCoeffs = computeGammaCoefficients(gamma,N);
  // Compute matrix entries
    for(int k=0; k<=N; k++){
      for(int l=0; l<=N; l++){
	// Compute double quadrature
	for(int qp1=0; qp1<4*N; qp1++){
	  auto s = TR_points(qp1);	  
	  for(int qp2=0; qp2<4*N; qp2++){
	    auto t = TR_points(qp2);
	    // compute argument of the log in hat{k}
	    auto evalCS = evaluateCosSinUn(t,s, N);
	    auto aux = (-(gammaCoeffs.block(0,0,N+1,2)).transpose()*evalCS.segment(0,N+1)
			    +(gammaCoeffs.block(N+1,0,N,2)).transpose()*evalCS.segment(N+1,N) );
	    // add contribution to the different matrices
	    Mcc(k,l)     += -TR_w/(4.*M_PI) * log(aux.squaredNorm())*cos(k*t)*cos(l*s);

	    if(k>0){
	      Msc(k-1,l) += -TR_w/(4.*M_PI) * log(aux.squaredNorm())*sin(k*t)*cos(l*s);
	    }

	    if(k>0 && l>0){
	      Mss(k-1,l-1) += -TR_w/(4.*M_PI) * log(aux.squaredNorm())*sin(k*t)*sin(l*s);
	    }	    	    
	  }// end loop over quadrature points
	}
      }// end for loop over l
    } // end for loop over k

    // Assemble big matrix
    Eigen::MatrixXd M(2*N+1, 2*N+1);
    std::cout << " Blocks " << std::endl;
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
      auto z = TR_points(qp);
      // First part (with $beta_k^c$)
      RHS(k) += TR_w * g(gamma(z)) * cos(k*z);
      // Second part (with $beta_k^s$)
      if(k>0){
	RHS(k+N) += TR_w * g(gamma(z)) * sin(k*z);
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
  // Build RHS
  std::cout << " Assemble RHS" << std::endl; 
  Eigen::VectorXd RHS = computeG(gamma, g, N);
  std::cout << " Solve " << std::endl;
  // Declare GMRES solver
  Eigen::ConjugateGradient< Eigen::MatrixXd > solver;
  solver.setMaxIterations( M.cols() );
  solver.setTolerance( 1.e-5 );
  // initialize the solver
  Eigen::MatrixXd LHS = AmM +M;
  solver.compute( LHS );
  Eigen::VectorXd sol = solver.solve(RHS);
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

int main() {

 std::cout << "==========  Test integration  ==========" << std::endl;
 int N1 = 30;
  Eigen::MatrixXd TR_points(N1, 2);  double TR_w;
  std::tie(TR_points,TR_w) = PeriodicTrapRule(N1);
  double Qint1 = 0.;  double Qintx = 0.;  double Qintcos = 0.;
  for(int qp=0; qp<N1; qp++){
    auto z = TR_points(qp);
    Qint1 += TR_w;
    Qintx += TR_w*z;
    Qintcos += TR_w*cos(z);
  }
  std::cout << " int 1 :" << Qint1 << " vs " << 2*M_PI << "\n int x :" << Qintx
	    << " vs " << 2*M_PI*M_PI << "\n int cosx :"  << Qintcos
	    << " vs 0" << std::endl;
  std::cout << "========================================" << std::endl << std::endl;

  
  std::cout << "============  Test Coefficients for S(t) = (cos(t), sin(t))  ============="
	    << std::endl;
  int N = 3;
    std::function<Eigen::Vector2d(const double&)> S = [](const double& t){
    Eigen::Vector2d res;
    res << cos(t) , sin(t);
    return res;
  };
  auto SCoeffs = computeGammaCoefficients(S,N);
  auto evalFun = evaluateCosSinUn(0.1, 0, N);
  Eigen::Vector2d diffS = (-SCoeffs.block(0,0,N+1,2)).transpose()*evalFun.segment(0,N+1)
              +(SCoeffs.block(N+1,0,N,2)).transpose()*evalFun.segment(N+1,N);
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
    res << cos(t) + 0.65*cos(2*t) - 0.65 , 1.5*sin(t);
    return res;
  };
  auto gammaCoeffs = computeGammaCoefficients(gamma,N);
  Eigen::Vector2d diffGamma = (-gammaCoeffs.block(0,0,N+1,2)).transpose()
                               *evalFun.segment(0,N+1)
                      +(gammaCoeffs.block(N+1,0,N,2)).transpose()*evalFun.segment(N+1,N);
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
  
  Eigen::VectorXd solC = solveBIE(S, gC, 5);
  double solCEval = reconstructRho(solC, M_PI, Sprime);
  std::cout << "eval at PI: " << solCEval
	    << " vs " << gC(S(M_PI))*2 << std::endl;
  std::cout << "====================================================================="
	    << std::endl << std::endl;

  // THE ISSUE IS IN M!!
  
  std::cout << "================  Test L2-norm  ================"
	    << std::endl;
    Eigen::VectorXd coeffToy(7);
  coeffToy << 0,1,0,0,0,0,0;
  double l2norm = L2norm(Sprime, coeffToy);
  std::cout << "L2 norm of cos(t) : " << l2norm << " vs " << std::sqrt(M_PI)
	    << std::endl
	    << "================================================"
	    << std::endl << std::endl;
    
    /*
  // SYSTEM
  std::function<Eigen::Vector2d(const double&)> gammaprime = [](const double& t){
    Eigen::Vector2d res;
    res << -sin(t) - 1.3*sin(2*t) , 1.5*cos(t);
    return res;
  };

  std::function<double(const Eigen::Vector2d&)> g = [](const Eigen::Vector2d& X){
    return sin(X(0))*sinh(X(1));
  };

  Eigen::VectorXd sol = solveBIE(gamma, g, 10);
  
  //std::cout << sol << std::endl;

  Eigen::Vector2d T({0.5, 0.2});
  double solEval = reconstructRho(sol, 0.1, gammaprime);
  std::cout << "eval at (0.5, 0.2) : " << solEval
	    << " vs " << g(gamma(0.1)) << std::endl;
  
    */  
    
  return 0;

}
