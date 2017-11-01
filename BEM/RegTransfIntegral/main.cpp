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
    v[0]=a;
    v[1]=b;
    v[2]=c;
  }

  Eigen::Vector3d getVertex(int i)const{
    assert(i>=0 && i<=2);
    return v[i]; 
  }
  
};


//-----------------------------------------------------------------------------
Eigen::Vector3d transform2BarycentricCoordinates(const Eigen::Vector3d& a,
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
  
  Eigen::Vector3d coeffs = transform2BarycentricCoordinates(a,b,c,x);
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
    // Use law of cosines to compute angle at b
    double C = b.norm();
    double B = c.norm();
    double A = (c-b).norm();
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
double distancePointToSegment(const Eigen::Vector3d& p, const Eigen::Vector3d& a,
			      const Eigen::Vector3d& b)
{
  Eigen::Vector3d dab = b - a;    // Distance between a and b. 
  Eigen::Vector3d r; r.setZero(); // Distance between [a,b] and p.

  if (dab.norm() <= 1e-12) { 
    r = p - a; // If a=b, return distance of p to a
  }
  else {
    // Consider the parametrization of the segment [a,b] : gamma(t) = a + t(b-a)
    // The projection p* of p on this line segment is given on
    double t = (p-a).dot(dab) / (dab.squaredNorm());
    // Now, return r = p-p* for the different cases:
    if (t <= 1e-12) {
      r = p-a; // if t<0, p*=a
    }
    else if ((1-t) <= 1e-12) {
      r = p-b; // if t>1, p*=b
    }
    else {
        r = p - a - t*dab; // t is in [0,1] (case when p* belongs to [a,b])
    }
  }

  return r.norm();
}


//-----------------------------------------------------------------------------
void transformCoordinates(const Eigen::Vector3d& a, const Eigen::Vector3d& b,
			  const Eigen::Vector3d& c, Eigen::Vector2d& v1,
			  Eigen::Vector2d& v2){
  Eigen::MatrixXd M(3,2);
  M << b-a, c-a;
  Eigen::HouseholderQR< Eigen::MatrixXd > qr(M);
  Eigen::MatrixXd R = qr.householderQ();
  v1 << (b-a).norm(), 0;
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
  const Eigen::Vector3d& a = T.getVertex(0);
  const Eigen::Vector3d& b = T.getVertex(1);
  const Eigen::Vector3d& c = T.getVertex(2);

  // if projection is inside the triangle T
  if(onTria){
    // Check whether it is close to one of the sides
    double dist_ab = distancePointToSegment(xp, a, b);
    double dist_bc = distancePointToSegment(xp, b, c);
    double dist_ca = distancePointToSegment(xp, c, a);
    if(dist_ab < 1e-12){
      std::cout << "xp close to ab" << std::endl;
      // Split triangle on 2
      { // T1: (xp , c , a)
	Eigen::Vector2d v1,v2;
	transformCoordinates(xp, c, a, v1, v2);
	res += integrateTiSing(v1, v2, (x-xp).norm(), n);
      }
      {	// T2: (xp , b , c)
	Eigen::Vector2d v1,v2;
	transformCoordinates(xp, b, c, v1, v2);
	res += integrateTiSing(v1, v2, (x-xp).norm(), n);
      }	
    }
    else if(dist_bc < 1e-12){
      std::cout << "xp close to bc" << std::endl;
      // Split triangle on 2
      {// T1: (xp , a , b)
	Eigen::Vector2d v1,v2;
	transformCoordinates(xp, a, b, v1, v2);
	res += integrateTiSing(v1, v2, (x-xp).norm(), n);
      }
      {	// T2: (xp , c , a)
	Eigen::Vector2d v1,v2;
	transformCoordinates(xp, c, a, v1, v2);
	res += integrateTiSing(v1, v2, (x-xp).norm(), n);
      }	
    }
    else if(dist_ca < 1e-12){
      std::cout << "xp close to ca" << std::endl;
      // Split triangle on 2
      { // T1: (xp , b , c)
	Eigen::Vector2d v1,v2;
	transformCoordinates(xp, b, c, v1, v2);
	res += integrateTiSing(v1, v2, (x-xp).norm(), n);
      }
      {// T2: (xp , a , b)
	Eigen::Vector2d v1,v2;
	transformCoordinates(xp, a, b, v1, v2);
	res += integrateTiSing(v1, v2, (x-xp).norm(), n);
      }
	
    }    
    else{
      // If not, split triangle on 3
      {// T1 : (xp, a, b)
	Eigen::Vector2d v1,v2;
	transformCoordinates(xp, a, b, v1, v2);
	res += integrateTiSing(v1, v2, (x-xp).norm(), n);
      }
      {// T2 : (xp, b, c)
	Eigen::Vector2d v1,v2;
	transformCoordinates(xp, b, c, v1, v2);
	res += integrateTiSing(v1, v2, (x-xp).norm(), n);
      }
      {// T3 : (xp, c, a)
      	Eigen::Vector2d v1,v2;
	transformCoordinates(xp, c, a, v1, v2);
	res += integrateTiSing(v1, v2, (x-xp).norm(), n);
      }	
    }
  }
  else{
    // still check that xp is not too close to T
    // if yes, treat it specially
    std::cout << " nau nau " << std::endl;
    // if not, we don't have a singularity and we can integrate
    // using gauss-legendre without any particular treatment.
    
  }
  
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
  double dist = distancePointToSegment(x, b1, c1);
  std::cout << "distance " << dist << std::endl;

  Eigen::Vector3d a, b, c;
  a << 1., 1., 1.;
  b << 2., 1., 0.;
  c << 0., 1., 0.;
  TriaPanel T(a, b, c);
  double It = integrateTSing(T, a, 6);
  
  return 0;

}
