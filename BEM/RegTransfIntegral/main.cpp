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


#if SOLUTION
//-----------------------------------------------------------------------------
/* SAM_LISTING_BEGIN_0 */
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
/* SAM_LISTING_END_0 */
#endif // SOLUTION


//-----------------------------------------------------------------------------
/* SAM_LISTING_BEGIN_1 */
bool ProjectOnTria(const TriaPanel& T, const Eigen::Vector3d& x,
		   Eigen::Vector3d& xp){
#if SOLUTION
  // Retrieve vertices from triangle T
  const Eigen::Vector3d& a = T.getVertex(0);
  const Eigen::Vector3d& b = T.getVertex(1);
  const Eigen::Vector3d& c = T.getVertex(2);
  
  Eigen::Vector3d coeffs = transformx2BarCoord(a,b,c,x);
  xp = coeffs(0)*a + coeffs(1)*b + coeffs(2)*c;

  // Check whether the point is inside of T or not
  if(coeffs(0)>0 && 1>coeffs(0) && coeffs(1)>0 && 1>coeffs(1)
     && coeffs(2)>0 && 1>coeffs(2))
    return true;
  else
    return false;
#else // TEMPLATE
  // TODO: Implement your code
  return true;
#endif // TEMPLATE
}
/* SAM_LISTING_END_1 */


//-----------------------------------------------------------------------------
#if SOLUTION
/* SAM_LISTING_BEGIN_2 */
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
    // Use law of cosines to compute angle at b and a (phimax)
    phib = std::acos( (A*A + C*C - B*B)/(2*A*C) );
    phiMax = std::acos( (B*B + C*C - A*A)/(2*B*C) );
    // Set radius when phi=0
    R0 = C;
    if( fabs(M_PI - phib-phiMax) < 1e-08 || fabs(M_PI + phib+phiMax) < 1e-08){
      std::cout << " Warning!! phib+phiMax too close to +-pi " << std::endl;
    }
    if( fabs(M_PI - phib) < 1e-08 || fabs(M_PI + phib) < 1e-08){
      std::cout << " Warning!! phib too close to pi " << std::endl;
    }
    if( fabs(phib) < 1e-08){
      std::cout << " Warning!! phib too close to 0 " << std::endl;
    }
  }

  double getPhiMax(){
    return phiMax;
  }

  double operator()(const double& phi){
    if( fabs(M_PI - phib-phi) < 1e-08 || fabs(M_PI + phib+phi) < 1e-08){
      std::cout << " Warning!! phib+phi too close to +-pi " << std::endl;
    }
    if( fabs(phib+phi) < 1e-08){
      std::cout << " Warning!! phib+phi too close to 0 " << std::endl;
    }
    // Apply formula derived from law of sines.
    return R0*sin(phib)/sin(phi+phib); 
  }

};
/* SAM_LISTING_END_2 */
#endif // SOLUTION


//-----------------------------------------------------------------------------
/* SAM_LISTING_BEGIN_3 */
double integrateTiSing(const Eigen::Vector2d& b, const Eigen::Vector2d& c,
		       double zeta, int n=6){
  double res=0;
#if SOLUTION
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
#else // TEMPLATE
    // TODO: Implement your code
#endif // TEMPLATE
  return res;
}
/* SAM_LISTING_END_3 */


#if SOLUTION
//-----------------------------------------------------------------------------
/* SAM_LISTING_BEGIN_4 */
void transformCoordinates(const Eigen::Vector3d& a, const Eigen::Vector3d& b,
			  const Eigen::Vector3d& c, Eigen::Vector2d& v1,
			  Eigen::Vector2d& v2){
  Eigen::MatrixXd M(3,2);
  M << b-a, c-a;
  // Use QR decomposition to get vectors for b and c on new basis
  Eigen::HouseholderQR< Eigen::MatrixXd > qr(M);
  Eigen::Matrix2d R = qr.matrixQR().template triangularView<Eigen::Upper>();
  // For simplicity, the rest of the code assumes that the transformation is on
  // the first quadrant. We check whether this is satisfied and flip the triangle
  // if not.
  if( R(0,0) < 0){    R.row(0) *= -1;  }
  if( R(1,1) < 0){    R.row(1) *= -1;  }
  // Finally, express vectors on new basis.
  v1 = R.col(0); // v1 for b
  v2 = R.col(1); // v2 for c
}
/* SAM_LISTING_END_4 */
#endif // SOLUTION


//-----------------------------------------------------------------------------
/* SAM_LISTING_BEGIN_5 */
double integrateTSing(const TriaPanel& T, const Eigen::Vector3d& x,
		      int n=6){
  double res = 0.;
#if SOLUTION
  // get projection of x on triangle's plane
  Eigen::Vector3d xp;
  bool onTria = ProjectOnTria(T, x, xp);
  // Retrieve vertices from triangle T
  std::array<Eigen::Vector3d, 3> vT = T.v;
  // Get barycentric coordinates of x
  const Eigen::Vector3d& coeff = transformx2BarCoord(vT[0],vT[1],vT[2],x);

  // We set the tolerance~sqrt(machine precision) to avoid triangles that might 
  // behave numerically as degenerate.
  double tol = 1e-08;
  // if projection is inside the triangle T (interior point)
  
  if(onTria){
   // Anyways check whether xp is close to the boundary of T
    bool boundary = false;
    // - Check if close to a vertex
    for(int i=0; i<3; i++){      
      if( fabs(1 - coeff(i)) < tol && !boundary){
	// Then we are close to vertex i. We don't split the triangle T but move the
	// given vertex to the origin and then integrate. Note that we can do it 
	// analytically when x=xp (see 1.11.c). For generality we use quadrature
	// for all cases.
	Eigen::Vector2d v1,v2;
	transformCoordinates(xp, vT[(i+1)%3], vT[(i+2)%3], v1, v2);
	res = integrateTiSing(v1, v2, (x-xp).norm(), n);
        boundary = true;
      }
    }    
    // - Check if close to an edge
    for(int i=0; i<3; i++){
      if ( fabs(coeff(i)) < tol && !boundary){
	// Then we are close to the edge opposite to the vertex i. We split the
	// triangle T into two smaller triangles
	/* 
	   vT(i+2)
	   |\              T1 : { xp, vT(i+2), vT(i) }
	   | \             T2 : { xp, vT(i+1), vT(i) }
	   |  \
	   |   \ xp         
	   |T1 "\   
	   |  "  \
	   | " T2 \
	   |"______\
  	   vT(i)    vT(i+1)
	*/
	for(int k=0; k<2; k++){
	  // Transform coordinates
	  Eigen::Vector2d v1,v2;
	  transformCoordinates(xp, vT[(i-k+2)%3], vT[i], v1, v2);
	  Rfunction Rtheta(v1,v2);
	  int nq = n;
	  if(n<100){
	  int m = (int)(Rtheta.getPhiMax());
	  nq = std::max(n, m*n);
	}	
	  // Add result from integrating over the current small triangle
	  res += integrateTiSing(v1, v2, (x-xp).norm(), nq);
	}
        boundary = true;
      }
    }
    // - If not close to boundary, we proceed as the standard interior point
    if(!boundary){
      // In this case, xp in an interior point of T and we split the triangle
      // into three smaller triangles.
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
      for(int k=0; k<3; k++){
	// Transform coordinates
	Eigen::Vector2d v1,v2;
	transformCoordinates(xp, vT[k], vT[(k+1)%3], v1, v2);
	Rfunction Rtheta(v1,v2);
	int nq = n;
	if(n<100){
	  int m = (int)(Rtheta.getPhiMax());
	  nq = std::max(n, m*5);
	}	  
	// Add result from integrating over the current small triangle
	res += integrateTiSing(v1, v2, (x-xp).norm(), nq);
      }
    }
     
  } // end if onTria
  // If xp is not on T
  else{
    // xp outside of T could mean that it is at the boundary of T or in the exterior.
    bool exterior = true;
    // Check whether xp is close to one of the vertices
    for(int i=0; i<3; i++){      
      if( fabs(1 - coeff(i)) < tol && fabs(coeff((i+1)%3)) < tol
	  && fabs(coeff((i+2)%3))< tol ){
	// Then we are on vertex i. We proceed as above.
	Eigen::Vector2d v1,v2;
	transformCoordinates(xp, vT[(i+1)%3], vT[(i+2)%3], v1, v2);
	res = integrateTiSing(v1, v2, (x-xp).norm(), n);
        exterior = false;
      }
    }// end for loop for vertices
    
    // Check whether xp is close to an edge. 
    for(int i=0; i<3; i++){
      if ( fabs(coeff(i)) < tol && exterior ){
	// Then we are on the edge opposite to the vertex i. We split the
	// triangle T into two smaller triangles as explained above.
	for(int k=0; k<2; k++){
	  // Transform coordinates
	  Eigen::Vector2d v1,v2;
	  transformCoordinates(xp, vT[(i-k+2)%3], vT[i], v1, v2);
	  Rfunction Rtheta(v1,v2);
	  int nq = n;
	  if(n<100){
	    int m = (int)(Rtheta.getPhiMax());
	    nq = std::max(n, m*n);
	  }	
	  // Add result from integrating over the current small triangle
	  res += integrateTiSing(v1, v2, (x-xp).norm(), nq);
	}
	exterior = false;
      }
    }// end for loop for edges

    // If xp is not too close to T, we don't have a singularity and we can 
    // integrate using gauss-legendre without any particular treatment. However,
    // as we learnt in the lecture (1.4.182), the proximity of a singularity
    // will be "felt" by Gaussin quadrature, so we choose the number of points
    // accordingly.
    if(exterior){
      // find closest point to x which is on T.
      Eigen::Vector3d coeff2;
      for(int i=0; i<3; i++){
	if(coeff(i)<=0.){ coeff2(i) = 0.;}
	else if(coeff(i)>=1.){ coeff2(i) = 1.;}
	else{ coeff2(i) = coeff(i);}
      }
      Eigen::Vector3d xpT = coeff2(0)*vT[0] + coeff2(1)*vT[1] + coeff2(2)*vT[2];
      // take distance to that point
      double distT_x = (x - xpT).norm();
      // choose number of points according to the distance to X
      int nq = 5*std::max(1, 1 - (int)log(distT_x) );
      nq = std::max(nq, n); // don't loose accuracy!
      // Get quadrature points and weights for [0, 1]
      // to use tensor-Gauss quadrature
      Eigen::RowVectorXd xq,wq;
      std::tie(xq,wq) =  gauleg(0, 1, nq);
      // Define transformation from reference triangle K to T
      Eigen::MatrixXd Fk(3,2);
      Fk << vT[1]-vT[0], vT[2]-vT[0];
      // Note we are mapping from 2D to 3D! (K to T)
      auto phiK = [&vT, &Fk](const Eigen::Vector2d& z){
	Eigen::Vector3d phiKz = Fk*z + vT[0];
	return phiKz;
      };
      // Integrate
      double Iq = 0.;
      for(int l=0; l<nq; l++){
	// Compute quadrature point on unit square
	Eigen::Vector2d qp; qp << xq(l), xq(l)*(1.-xq(l));
	// Multiply quadrature weight by determinants coming from transformations
	double dy = wq(l)*(1.-xq(l))*sqrt(fabs((Fk.transpose()*Fk).determinant()));
	// add contribution
	Iq += 1./(x - phiK(qp)).norm()*dy;
      }
      res = Iq;
    } // end if not too close
    
  } // end if outside
#else // TEMPLATE
  // TODO: Implement your code
#endif // TEMPLATE
  return res;
}
/* SAM_LISTING_END_5 */

//------------------------------------------------------------------------------
void testIntegrateTSing(const int& n){
  // Test on reference triangle for x on vertices and edges (for which we 
  // computed the exact solution using Maple).
  Eigen::Vector3d ah, bh, ch;
  ah << 0., 0., 0.;
  bh << 1., 0., 0.;
  ch << 0., 1., 0.; 
  TriaPanel T0(ah, bh, ch);
  
  std::cout << " Testing for n = " << n << std::endl;

  std::cout  << " Error on a :"
	     << fabs(integrateTSing(T0, ah, n)  + sqrt(2)*log(sqrt(2)-1) )
	     << std::endl;

  std::cout << " Error on b :"
	    << fabs(integrateTSing(T0, bh, n) + log(sqrt(2)-1) )
	    << std::endl;

  std::cout << " Error on c :"
	    << fabs(integrateTSing(T0, ch, n) + log(sqrt(2)-1) )
	    << std::endl;

  std::cout << " Error on 0.5*(a+b) :" ;
  double exres = 1.676348272;
  std::cout << fabs(integrateTSing(T0, 0.5*(ah + bh), n) - exres )
	    << std::endl;
  
  std::cout << " Error on 0.5*(b+c) :" ;
  double exres2 = 0.5*( log(sqrt(2)+1) - 3*log(sqrt(2)-1) );
  std::cout << fabs(integrateTSing(T0, 0.5*(bh + ch), n) - exres2 )
	    << std::endl;
  
  std::cout << " Error on 0.5*(a+c) :" ;
  std::cout << fabs(integrateTSing(T0, 0.5*(ah + ch), n) - exres )
	    << std::endl  << std::endl;
  
}


//-----------------------------------------------------------------------------
int main() {

  // 1. Test convergence for integration over Ti
  Eigen::Vector2d b0, c0; b0<< 1., 0.; c0 << 0.5, 1.;
  double Iref = integrateTiSing(b0, c0, 0., 1000);
  Eigen::VectorXd error(19);
  Eigen::VectorXi N = Eigen::VectorXi::LinSpaced(19,1,19);
  for(int k=1; k<20; k++){
    double I = integrateTiSing(b0, c0, 0., k);
    error(k-1) = fabs(I - Iref);
  }
  // Output for plot
  std::ofstream out_errorQ("integrateTiSing_errors.txt");
  out_errorQ << std::setprecision(18) << error; 
  out_errorQ.close( );
  std::ofstream out_N("integrateTiSing_N.txt");
  out_N << N; 
  out_N.close( );

  // 2. Test integration routine on reference triangle (help for debugging)
  testIntegrateTSing(5);
  testIntegrateTSing(10);
  testIntegrateTSing(20);
  
  // 3. Test for given triangle moving the point x
  Eigen::Vector3d a, b, c;  
  a << 1., 1., 1.;
  b << 2., 1., 0.;
  c << 0., 1., 0.; 
  TriaPanel T(a, b, c);

  // We choose number of points to discretize [0, 2]
  int Nts = 101;
  Eigen::VectorXd tau = Eigen::VectorXd::LinSpaced(Nts, 0, 2);
  // Create vectors to save the results
  Eigen::VectorXd intx5(Nts), intx10(Nts), intx20(Nts), intx1000(Nts);
  for(int ts=0; ts<Nts; ts++){
    Eigen::Vector3d xtau = a + tau(ts)*( 0.5*(b+c) - a);
    std::cout << " tau " << tau(ts) << " is ";
    intx5(ts) = integrateTSing(T, xtau, 5);
    intx10(ts) = integrateTSing(T, xtau, 10);
    intx20(ts) = integrateTSing(T, xtau, 20);
    intx1000(ts) = integrateTSing(T, xtau, 1000);
    
  }
  
  // Output errors for plot
  std::ofstream out_Intx5("intx5.txt");
  out_Intx5 << std::setprecision(18) << (intx1000 - intx5).array().abs(); 
  out_Intx5.close( );
  std::ofstream out_Intx10("intx10.txt");
  out_Intx10 << std::setprecision(18) << (intx1000 - intx10).array().abs();  
  out_Intx10.close( );
  std::ofstream out_Intx20("intx20.txt");
  out_Intx20 << std::setprecision(18) << (intx1000 - intx20).array().abs(); 
  out_Intx20.close( );
  std::ofstream out_dt("dt.txt");
  out_dt << std::setprecision(18) << tau; 
  out_dt.close( );
  std::ofstream out_IntxOK("intx1000.txt");
  out_IntxOK << std::setprecision(18) << intx1000 ;
  out_IntxOK.close( );
  
    
  return 0;

}
