//// 
//// Copyright (C) 2017 SAM (D-MATH) @ ETH Zurich
//// Author(s): curzuato < > 
//// Contributors:  dcasati 
//// This file is part of the AdvNumCSE repository.
////
#include <iostream>
#include <fstream>
#include <istream>
#include <iomanip>
#include <cmath>
#include <Eigen/Dense>
#include "gauleg.hpp"


//-----------------------------------------------------------------------------
struct TriaPanel{
  // Array of 3d-vectors containing the triangle's vertices
  std::array<Eigen::Vector3d, 3> v;
  TriaPanel( const Eigen::Vector3d& a, const Eigen::Vector3d& b,
	     const Eigen::Vector3d& c){
    // Initialize array of vertices with the corresponding points
    v[0]=a;  v[1]=b;  v[2]=c;
  }

  Eigen::Vector3d getVertex(int i)const{
    assert(i>=0 && i<=2);
    return v[i]; 
  }
  
};




//-----------------------------------------------------------------------------
bool ProjectOnTria(const TriaPanel& T, const Eigen::Vector3d& x,
		   Eigen::Vector3d& xp){
  // TODO: Implement your code
  return true;
}


//-----------------------------------------------------------------------------


//-----------------------------------------------------------------------------
double integrateTiSing(const Eigen::Vector2d& b, const Eigen::Vector2d& c,
		       double zeta, int n=6){
  double res=0;
    // TODO: Implement your code
  return res;
}




//-----------------------------------------------------------------------------
double integrateTSing(const TriaPanel& T, const Eigen::Vector3d& x,
		      int n=6){
  double res = 0.;
  // TODO: Implement your code
  return res;
}



//-----------------------------------------------------------------------------
int main() {
  
  Eigen::Vector2d b0, c0; b0<< 1., 0.; c0 << 0.5, 1.;
  double Iref = integrateTiSing(b0, c0, 0., 1000);
  Eigen::VectorXd error(19);
  Eigen::VectorXi N = Eigen::VectorXi::LinSpaced(19,1,19);
  for(int k=1; k<20; k++){
    double I = integrateTiSing(b0, c0, 0., k);
    error(k-1) = fabs(I - Iref);
  }
  //output for plot
  std::ofstream out_errorQ("integrateTiSing_errors.txt");
  out_errorQ << std::setprecision(18) << error; 
  out_errorQ.close( );
  std::ofstream out_N("integrateTiSing_N.txt");
  out_N << N; 
  out_N.close( );

  Eigen::Vector3d b1, c1; b1<< 1., 0., 0.; c1 << 0.5, 1., 0.;
  Eigen::Vector3d x; x<< 0.75, 0.65, 0.;
  //double dist = distancePointToSegment(x, b1, c1);
  //std::cout << "distance " << dist << std::endl;
  /*
  // Test quadrature
  Eigen::Vector3d a, b, c;  
  a << 1., 1., 1.;
  b << 2., 1., 0.;
  c << 0., 1., 0.; 
  TriaPanel T(a, b, c);  
  int Nts = 100;
  Eigen::VectorXd tau = Eigen::VectorXd::LinSpaced(Nts, 0, 2);
  Eigen::VectorXd intx5(Nts), intx10(Nts), intx20(Nts);
  for(int ts=0; ts<Nts; ts++){
    Eigen::Vector3d xtau = a + tau(ts)*( 0.5*(b+c) - a);
    intx5(ts) = integrateTSing(T, xtau, 5);
    intx10(ts) = integrateTSing(T, xtau, 10);
    intx20(ts) = integrateTSing(T, xtau, 20);   
  }
  std::ofstream out_Intx5("intx5.txt");
  out_Intx5 << std::setprecision(18) << intx5; 
  out_Intx5.close( );
  std::ofstream out_Intx10("intx10.txt");
  out_Intx10 << std::setprecision(18) << intx10; 
  out_Intx10.close( );
  std::ofstream out_Intx20("intx20.txt");
  out_Intx20 << std::setprecision(18) << intx20; 
  out_Intx20.close( );
  std::ofstream out_dt("dt.txt");
  out_dt << std::setprecision(18) << tau; 
  out_dt.close( );
  */

  // Second test
  Eigen::Vector3d ah, bh, ch;  
  ah << 0., 0., 0.;
  bh << 1., 0., 0.;
  ch << 0., 1., 0.; 
  TriaPanel T0(ah, bh, ch);
  std::cout << " On a " << std::endl;
  std::cout << " Error for n = 5 "
    << fabs(integrateTSing(T0, ah, 5)  + sqrt(2)*log(sqrt(2)-1) )
	    << "\n Error for n = 10 "
    << fabs(integrateTSing(T0, ah, 10) + sqrt(2)*log(sqrt(2)-1)  )
	    << "\n Error for n = 20 "
    << fabs(integrateTSing(T0, ah, 20) + sqrt(2)*log(sqrt(2)-1)  )
	    << std::endl;

  std::cout << " On b " << std::endl;
  std::cout << " Error for n = 5 "    << fabs(integrateTSing(T0, bh, 5) + log(sqrt(2)-1) )
	    << "\n Error for n = 10 " << fabs(integrateTSing(T0, bh, 10) + log(sqrt(2)-1) )
	    << "\n Error for n = 20 " << fabs(integrateTSing(T0, bh, 20) + log(sqrt(2)-1) )
	    << std::endl;

  std::cout << " On c " << std::endl;
  std::cout << " Error for n = 5 "    << fabs(integrateTSing(T0, ch, 5) + log(sqrt(2)-1) )
	    << "\n Error for n = 10 " << fabs(integrateTSing(T0, ch, 10) + log(sqrt(2)-1) )
	    << "\n Error for n = 20 " << fabs(integrateTSing(T0, ch, 20) + log(sqrt(2)-1) )
	    << std::endl;

  std::cout << " On 0.5*(a+c) " << std::endl;
  double exres = 1.676348272;
  std::cout << " Error for n = 5 "
	    << fabs(integrateTSing(T0, 0.5*(ah + ch), 5) - exres )
	    << "\n Error for n = 10 "
    	    << fabs(integrateTSing(T0, 0.5*(ah + ch), 10) - exres )
	    << "\n Error for n = 20 "
    	    << fabs(integrateTSing(T0, 0.5*(ah + ch), 20) - exres )
	    << std::endl;
  
  return 0;

}
