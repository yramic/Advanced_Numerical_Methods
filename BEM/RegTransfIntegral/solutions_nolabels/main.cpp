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
Eigen::Vector3d transformx2BarCoord(const Eigen::Vector3d& a,
				    const Eigen::Vector3d& b,
				    const Eigen::Vector3d& c,
				    const Eigen::Vector3d& x){
  Eigen::Vector3d coeffs;
  // Define vectors to compute projection on barycentric coordinates
  Eigen::Vector3d u = b-a;
  Eigen::Vector3d v = c-a;
  Eigen::Vector3d w = x-a;
  Eigen::Vector3d n = u.cross(v);  
  // Compute barycentric coordinates of xp in the plane of the triangle T
  coeffs(1) = (w.cross(v)).dot(n)/n.squaredNorm();
  coeffs(2) = (u.cross(w)).dot(n)/n.squaredNorm();
  coeffs(0) = 1 - coeffs(1) - coeffs(2);
  return coeffs;
}


//-----------------------------------------------------------------------------
bool ProjectOnTria(const TriaPanel& T, const Eigen::Vector3d& x,
		   Eigen::Vector3d& xp){
  // Retrieve vertices from triangle T
  const Eigen::Vector3d& a = T.getVertex(0);
  const Eigen::Vector3d& b = T.getVertex(1);
  const Eigen::Vector3d& c = T.getVertex(2);
  
  Eigen::Vector3d coeffs = transformx2BarCoord(a,b,c,x);
  xp = coeffs(0)*a + coeffs(1)*b + coeffs(2)*c;

  // Check whether the point is inside of T or not
  if(coeffs(0)>=0 && 1>=coeffs(0) && coeffs(1)>=0 && 1>=coeffs(1)
     && coeffs(2)>=0 && 1>=coeffs(2))
    return true;
  else
    return false;
}


//-----------------------------------------------------------------------------
struct Rfunction{
  double phib;
  double R0;
  double phiMax;

  Rfunction(const Eigen::Vector2d& b, const Eigen::Vector2d& c){
    /* 
            c 
            /\
           /  \
          /    \         (notation for vertices and 
      B  /      \ A        sides in this code)
        /        \
       /          \
      /____________\
     a              b
            C
    */
    double C = b.norm();
    double B = c.norm();
    double A = (c-b).norm();
    // Use law of cosines to compute angle at b
    phib = std::acos( (A*A + C*C - B*B)/(2*A*C) );
    // Use law of sines to compute maximum angle (angle at a)
    phiMax = std::asin( A*sin(phib)/B);
    // Set radius when phi=0
    R0 = C;
  }

  double getPhiMax(){
    return phiMax;
  }

  double operator()(const double& phi){
    // Apply formula derived from law of sines.
    return R0*sin(phib)/sin(phi+phib); 
  }

};


//-----------------------------------------------------------------------------
double integrateTiSing(const Eigen::Vector2d& b, const Eigen::Vector2d& c,
		       double zeta, int n=6){
  double res=0;
  // Initialize R function for current triangle;
  Rfunction Rtheta(b,c);
  // Retrieve maximum angle
  double phimax = Rtheta.getPhiMax();
  // Get quadrature points and weights for [0, phimax]
  Eigen::RowVectorXd xq,wq;
  std::tie(xq,wq) =  gauleg(0, phimax, n);
  // Integrate
  for(int k=0; k<n; k++){
    res += (std::sqrt(Rtheta(xq(k))*Rtheta(xq(k))+zeta*zeta) - zeta)*wq(k);
  }
  return res;
}


//-----------------------------------------------------------------------------
void transformCoordinates(const Eigen::Vector3d& a, const Eigen::Vector3d& b,
			  const Eigen::Vector3d& c, Eigen::Vector2d& v1,
			  Eigen::Vector2d& v2){
  std::cout << " currently transforming : " << a.transpose() << " , "
	    << b.transpose() << " , " << c.transpose() << std::endl;
  Eigen::MatrixXd M(3,2);
  M << b-a, c-a;
  Eigen::HouseholderQR< Eigen::MatrixXd > qr(M);
  Eigen::MatrixXd R = qr.matrixQR().template triangularView<Eigen::Upper>();
  v1 = R.col(0);
  v2 = R.col(1);
}


//-----------------------------------------------------------------------------
double integrateTSing(const TriaPanel& T, const Eigen::Vector3d& x,
		      int n=6){
  double res = 0.;
  // get projection of x on triangle's plane
  Eigen::Vector3d xp;
  bool onTria = ProjectOnTria(T, x, xp);
  // Retrieve vertices from triangle T
  std::array<Eigen::Vector3d, 3> vT = T.v;
  // Get barycentric coordinates of x
  const Eigen::Vector3d& coeff = transformx2BarCoord(vT[0],vT[1],vT[2],x);
  std::cout << " I am " << coeff.transpose() ;

  // if projection is inside the triangle T
  if(onTria){
    bool interior = true;
    // Check whether xp is close to one of the vertices
    for(int i=0; i<3; i++){      
      if( fabs(1 - coeff(i)) < 1e-12 ){
	// Then we are on vertex i. We don't split the triangle T but move the
	// given vertex to the origin and then integrate. Note that we can do it 
	// analytically when x=xp (see 1.11.c). For generality we use quadrature
	// for all cases.
	std::cout << ". I am at vertex " << i << std::endl;
	Eigen::Vector2d v1,v2;
	transformCoordinates(xp, vT[(i+1)%3], vT[(i+2)%3], v1, v2);
	std::cout << v1.transpose() << " , " << v2.transpose() << std::endl;
	res = integrateTiSing(v1, v2, (x-xp).norm(), n);
	if(!interior){
	  std::cout << " weird " << i << std::endl;
	}
	interior = false;
      }
    }// end for loop for vertices
    
    // Check whether xp is close to an edge
    for(int i=0; i<3; i++){
      if ( fabs(coeff(i)) < 1e-12 && interior ){
	// Then we are on the edge opposite to the vertex i. We split the
	// triangle T into two smaller triangles
	/* 
	   vT(i+2)
	     |\              T1 : { xp, vT(i+2), vT(i+1) }
	     | \             T2 : { xp, vT(i)  , vT(i+2) }
	     |  \
             |   \ xp         
	     |T1 "\   
	     |  "  \
	     | " T2 \
	     |"______\
  	   vT(i)    vT(i+1)
	*/
	std::cout << ". I am at edge opposite to vertex " << i << std::endl;
	for(int k=0; k<2; k++){
	  // Transform coordinates
	  Eigen::Vector2d v1,v2;
	  transformCoordinates(xp, vT[(i+k+2)%3], vT[(i+k)%3], v1, v2);
	  std::cout << v1.transpose() << " , " << v2.transpose() << std::endl;
	  // Add result from integrating over the current small triangle
	  res += integrateTiSing(v1, v2, (x-xp).norm(), n);
	}
	if(!interior){
	  std::cout << " weirdo " << i << std::endl;
	}
	interior = false;
      }
    }// end for loop for edges
    
    if(interior){ 
      // In this case, xp in an interior point of T and we split the triangle
      // into three smaller triangles
      /* 
	          vT(2)
	           /"\
                  / " \
	         /  "  \          T1 : { xp, vT(0), vT(1) }
	        /  xp   \         T2 : { xp, vT(1), vT(2) }
	       /  "  "   \        T3 : { xp, vT(2), vT(0) }
	      / "      "  \
	     /"___________"\
	   vT(0)          vT(1) 
      */
      std::cout << ", I am interior " << std::endl;
      for(int k=0; k<3; k++){
	// Transform coordinates
	Eigen::Vector2d v1,v2;
	transformCoordinates(xp, vT[k], vT[(k+1)%3], v1, v2);
	// Add result from integrating over the current small triangle
	res += integrateTiSing(v1, v2, (x-xp).norm(), n);
      }
    } 
  } // end if onTria
  // If xp is not on T
  else{
    std::cout << ". I am outside of T ";
    bool close2T = false;
    /*
    // Still check that xp is not too close to T (for this we do as before)
    // Check whether xp is close to one of the vertices
    for(int i=0; i<3; i++){
      if( fabs(1 - coeff(i)) < 1e-12 ){
	Eigen::Vector2d v1,v2;
	transformCoordinates(xp, vT[(i+1)%3], vT[(i+2)%3], v1, v2);
	res = integrateTiSing(v1, v2, (x-xp).norm(), n);
        close2T = true;
	std::cout << "wii" << std::endl;
      }
    }
    // Check whether xp is close to an edge
    for(int i=0; i<3; i++){
      if ( fabs(coeff(i)) < 1e-12 ){
	for(int k=0; k<2; k++){
	  Eigen::Vector2d v1,v2;
	  transformCoordinates(xp, vT[(i+k+2)%3], vT[(i+k)%3], v1, v2);
	  res += integrateTiSing(v1, v2, (x-xp).norm(), n);
	}
	std::cout << "wii?" << std::endl;
        close2T = true;
      }
    }
    */

    // If xp is not too close from T, we don't have a singularity and we can 
    // integrate using gauss-legendre without any particular treatment. However,
    // as we learnt in the lecture (1.4.182), the proximity of a singularity
    // will be "felt" by Gaussina quadrature, so we choose the number of points
    // accordingly.
    if(!close2T){
      // find closest point to x which is on T.
      Eigen::Vector3d coeff2;
      for(int i=0; i<3; i++){
	if(coeff(i)<=0)
	  coeff2(i) = 0;
	if(coeff(i) >=1)
	  coeff2(i) = 1;
	else
	  coeff2(i) = coeff(i);
      }
      Eigen::Vector3d xpT = coeff2(0)*vT[0] + coeff2(1)*vT[1] + coeff2(2)*vT[2];
      // take distance to that point
      double distT_x = (x - xpT).norm();
      // choose number of points accordingly
      int nq = n* std::max(1., 1+log(2*distT_x));
      std::cout << ". Using " << nq << " quadrature points" << std::endl;
      // Get quadrature points and weights for [0, 1]
      Eigen::RowVectorXd xq,wq;
      std::tie(xq,wq) =  gauleg(0, 1, nq);
      // Transform T to reference triangle K and then to unit square in order
      Eigen::MatrixXd Fk(3,2);
      Fk << vT[1]-vT[0], vT[2]-vT[0];
      // Note we are mapping from 2D to 3D! (K to T)
      auto phiK = [&vT, &Fk](const Eigen::Vector2d& z){
	Eigen::Vector3d phiKz = Fk*z + vT[0];
	return phiKz;
      };
      // to use tensor-Gauss quadrature
      double Iq = 0.;
      for(int l=0; l<nq; l++){
	// Compute quadrature point on unit square
	Eigen::Vector2d qp;
	qp << xq(l), xq(l)*(1.-xq(l));
	// Multiply quadrature weight by determinants coming from transformations
	double dy = wq(l)*(1.-xq(l))*sqrt((Fk.transpose()*Fk).determinant());
	// add contribution
	Iq += 1./(x - phiK(qp)).norm()*dy;
      }
      res = Iq;
    } // end if not too close
    
  } // end if outside
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
