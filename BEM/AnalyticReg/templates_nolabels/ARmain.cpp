//// 
//// Copyright (C) 2017 SAM (D-MATH) @ ETH Zurich
//// Author(s): dcasati <daniele.casati@sam.math.ethz.ch> 
//// Contributors: curzuato
//// This file is part of the AdvNumCSE repository.
////
#include <iostream>
#include <fstream>
#include <istream> 
#include <cmath>
#include <Eigen/Dense>
// GMRES
#include <unsupported/Eigen/IterativeSolvers>



Eigen::MatrixXd computeAminusM(int N){
  Eigen::DiagonalMatrix<double, Eigen::Dynamic> AmM(2*N+1);

  // TODO: Compute A-M
  return AmM;
}



template <typename PARAM>
Eigen::MatrixXd computeM(const PARAM& gamma, int N){

  // Assemble big matrix
  Eigen::MatrixXd M(2*N+1, 2*N+1);

    // TODO: Compute M using quadrature

  return M;
}


template <typename PARAM, typename FUNC>
Eigen::VectorXd computeG(const PARAM& gamma, const FUNC& g, int N){
  // Initialize right hand side vector
  Eigen::VectorXd RHS(2*N+1);  RHS.setZero();
    // TODO: Compute RHS
  return RHS;
}


template <typename PARAM, typename FUNC>
Eigen::VectorXd solveBIE(const PARAM& gamma, const FUNC& g, int N){
    // TODO: Build BIE system and solve it
  Eigen::VectorXd sol(2*N+1);
  
  return sol;
}

 
template <typename PARAMDER>
double L2norm(const PARAMDER& gammaprime, const Eigen::VectorXd& coeffs){
  double res=0.;
    // TODO: reconstruct the function and compute its L2norm
  
  return std::sqrt(res);
}



int main() {

  int N = 9;
  //----------------------------------------------------------------------------
  std::cout << "===========  Test integration  ===========" << std::endl;
  double Qint1 = 0., Qintcos = 0., Qintlogcos2 = 0.;
  // TODO: You may test your quadrature by integrating 1, cos(t) and log(cos(t)^2 + 1)
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
  std::function<Eigen::Vector2d(const double&)> Sprime = [](const double& t){
    Eigen::Vector2d res;
    res << -sin(t),  cos(t);
    return res;
  };

  
  // TODO: You may test your computation of the coefficients for S(t) and the
  // difference (S(0.1)-S(0))/S(0.1)-S(0).
  Eigen::Vector2d diffS; diffS.setZero();
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
  // TODO: You may test your computation of the coefficients for S(t) and the
  // difference (S(0.1)-S(0))/S(0.1)-S(0).
  Eigen::Vector2d diffGamma; diffGamma.setZero();
  Eigen::Vector2d exDiffgamma = gamma(0.1) - gamma(0);
  std::cout << "(gamma (0.1) - gamma(0))//S(0.1)-S(0) = " << diffGamma.transpose()
	    << " vs " << exDiffgamma.transpose()/exDiffS.norm() << std::endl;
  std::cout << "============================================================="
	    << "=============" << std::endl << std::endl;
  


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


  Eigen::VectorXi Nall(9); Nall<< 3,5,7,9,11,13,15,17,19;
  Eigen::VectorXd error(9); error.setZero();
  Eigen::Vector2d T({0.5,0.2});
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
  
  std::cout << "DISCLAIMER : This code is still not working! " << std::endl;
    
  return 0;

}
