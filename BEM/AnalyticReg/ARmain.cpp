#include <iostream>
#include <cmath>
#include <Eigen/Dense>
// GMRES
#include <unsupported/Eigen/IterativeSolvers>
#include "PeriodicTrapezoidalQR.hpp"


Eigen::MatrixXd computeAminusM(int N){
  Eigen::DiagonalMatrix<double, Eigen::Dynamic> AmM(2*N+1);
  Eigen::VectorXd aux(2*N+1);
  aux << 0, Eigen::VectorXd::LinSpaced(N,1,N).cwiseInverse()*(M_PI/2.),
    Eigen::VectorXd::LinSpaced(N,1,N).cwiseInverse()*(M_PI/2.);
  AmM = aux.asDiagonal();
  return AmM;
}

template <typename PARAM>
Eigen::MatrixXd computeGammaCoefficients(const PARAM& gamma, int N){
  // Get quadrature points and weight
  Eigen::MatrixXd TR_points(N, 2);
  double TR_w;
  std::tie(TR_points,TR_w) = PeriodicTrapRule(N);
  
  Eigen::MatrixXd coeffs(N+1,2);
  coeffs.setZero();

  for(int k=0; k<=N; N++){
    // iterate over quadrature points
    for(int qp=0; qp<N; qp++){
      //auto z = TR_points.row(qp);
      auto z = TR_points(qp);
    // compute coeficcients ak
      coeffs(k,0) += 1/M_PI*TR_w*gamma(z)*cos(k*z);

      if(k>0){
	// compute coefficients bk
	coeffs(k,1) += 1/M_PI*TR_w*gamma(z)*sin(k*z);
      }
    }// end for qp

  }// end for coeffs

  return coeffs;  
}

Eigen::MatrixXd evaluateCosSin(double z, int N){
  Eigen::MatrixXd res(N+1,2);
  res.setZero();
  for(int k=0; k<=N; N++){
    res(k,0) = cos(k*z);
    res(k,1) = sin(k*z);
  }

}

template <typename PARAM>
Eigen::MatrixXd computeM(const PARAM& gamma, int N){

  // Get quadrature points and weight
  Eigen::VectorXd TR_points(2*N);
  double TR_w;
  std::tie(TR_points,TR_w) = PeriodicTrapRule(2*N);

  // For readibility, build each separate block and then assemble the big matrix.
  Eigen::MatrixXd Mcc(N+1,N+1), Mss(N,N), Msc(N,N+1);
  // Initialize them to zero;
  Mcc.setZero(); Mss.setZero(); Msc.setZero(); 

  // In order to have a smooth integrand, we use the formula derived on 1.8.b
  // Get coefficients for Gamma in (1.0.42)
  auto gammaCoeffs = computeGammaCoefficients(gamma,N);

  // Compute matrix entries
    for(int k=0; k<=N; k++){
      for(int l=0; l<=N; l++){
	// Compute double quadrature
	for(int qp1=0; qp1<2*N; qp1++){
	  auto s = TR_points(qp1);
	  
	  for(int qp2=0; qp2<2*N; qp2++){
	    auto t = TR_points(qp2);

	    // compute argument of the log in \hat{k}
	    auto evalCS = evaluateCosSin(s-t, N);
	    auto aux = 0.5*(-evalCS.col(0).dot(gammaCoeffs.col(0))
			    + evalCS.col(1).dot(gammaCoeffs.col(1)) );

	    // add contribution to the different matrices
	    Mcc(k,l)     += -TR_w/(4.*M_PI) * log(aux*aux)*cos(k*t)*cos(l*s);

	    if(k>0){
	      Msc(k-0,l) += -TR_w/(4.*M_PI) * log(aux*aux)*sin(k*t)*cos(l*s);
	    }

	    if(k>0 && l>0){
	      Mss(k-0,l-0) += -TR_w/(4.*M_PI) * log(aux*aux)*sin(k*t)*sin(l*s); 
	    }
	    	    
	  }// end loop over quadrature points
	}
      }// end for loop over l
    } // end for loop over k

    // Assemble big matrix
    Eigen::MatrixXd M(2*N+1, 2*N+1);
    M.block(0  , 0  , N+1, N+1) = Mcc;
    M.block(0  , N+1, N+1, N  ) = Msc.transpose();
    M.block(N+1, 0  , N  , N+1) = Msc;
    M.block(N+1, N+1, N  , N  ) = Mss;

    return M;
}


template <typename PARAM>
double L2norm(const PARAM& gamma, const Eigen::VectorXd& coeffs){



}

template <typename PARAM, typename FUNC>
Eigen::VectorXd computeG(const PARAM& gamma, const FUNC& g, int N){

  // Initialize right hand side vector
  Eigen::VectorXd RHS(2*N+1);
  RHS.setZero();
  
  // Get quadrature points and weight
  Eigen::VectorXd TR_points(2*N);
  double TR_w;
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
	RHS(k+N+1) += TR_w * g(gamma(z)) * sin(k*z);
      }	
    }// end iteration over quadrature points
  } 
  
  return RHS;
}

template <typename PARAM, typename FUNC>
Eigen::VectorXd solveBIE(const PARAM& gamma, const FUNC& g, int N){
  // In order to compute A we do (A-M)+M
  Eigen::MatrixXd AmM = computeAminusM(N);
  Eigen::MatrixXd M   = computeM(gamma, N);

  // Build RHS
  Eigen::VectorXd RHS = computeG(gamma, g, N); 

  // Solve -> which solver do we use?
  Eigen::ConjugateGradient< Eigen::MatrixXd > solver;
  solver.setMaxIterations( M.cols() );
  solver.setTolerance( 1.e-5 );
  // initialize the solver
  Eigen::MatrixXd LHS = AmM+M ;
  solver.compute( LHS );
  Eigen::VectorXd sol = solver.solve(RHS); 

}

int main() {

  std::cout << " wii " << std::endl;
  

  return 0;

}
