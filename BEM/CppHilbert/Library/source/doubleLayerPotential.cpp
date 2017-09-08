///////////////////////////////////////////////////////////////////////////////
/// \file doubleLayerPotential.cpp
/// \brief This file provides functions to compute the Galerkin Matrix K given
///  by
///  \f[ K_{ij} = -\frac{1}{2\pi}\int_{Ei} \int_{supp \phi_j}\frac{<y-x,n>}{
///                |y-x|^2} \phi_j(y)ds_y ds_x \f]
///
///  This file contains only the implementation. For detailed documentation see
///  doubleLayerPotential.hpp
///
///  This file is part of the HILBERT program package for the numerical
///  solution of the Laplace equation with mixed boundary conditions by use of
///  BEM in 2D.
///
///  C++ adaptation for ANCSE17 of HILBERT V3.1 TUWien 2009-2013
////////////////////////////////////////////////////////////////////////////////
#include <cmath>
#include <cassert>

#include "constants.hpp"
#include "doubleLayerPotential.hpp"

extern "C" {
#include "gaussQuadrature.h"
}


//-----------------------------------------------------------------------------
double dlp(int k, const Eigen::Vector2d& p, const Eigen::Vector2d& q)
{
  double a = p.squaredNorm();  /* a = <p,p> */
  double b = 2 * p.dot(q);     /* b = 2 <p,q> */
  double c = q.squaredNorm();  /* c = <q,q> */
  double D = 4*a*c-b*b;
  double root_D = 0;
  double G0 = 0;
  double G1 = 0;

  assert(D>=-EPS*4*a*c);
  if (D > EPS*4*a*c)
    root_D = sqrt(D);
  else
    D = 0.0;

  if (D == 0.0)
  {
    G0 = 2./(c-a);
  }   
  else
  {
    if (fabs(c-a) < EPS*fabs(c))
      G0 = M_PI/root_D;
    else if (a < c)
      G0 = 2.*atan(root_D/(c-a))/root_D;
    else
      G0 = 2.*(atan(root_D/(c-a))+M_PI)/root_D;
  }

  if (k >= 1)
  {
    G1 = -b*G0;
    if (a+b+c > EPS*a)
      G1 += log(a+b+c);
      
    if (a-b+c > EPS*a)
      G1 -= log(a-b+c);

    G1 /= (2.*a);
    
    if (k == 2)
      return (2.-b*G1-c*G0)/a;
    
    return G1;
  }

  return G0;
}


//-----------------------------------------------------------------------------
void computeKij(double* I0, double* I1, double eta,
                const Eigen::Vector2d& a, const Eigen::Vector2d& b,
		const Eigen::Vector2d& c, const Eigen::Vector2d& d)
{
  int swap = 0;
  double hi = (b-a).squaredNorm(); /* hi = norm(b-a)^2 */
  double hj = (d-c).squaredNorm(); /* hj = norm(d-c)^2 */
  
  if ((hi-hj)/hi > EPS)
  {
    swap = 1;     /* swap Ej and Ei */
    std::swap(hi,hj);
  }

  if (eta == 0)   /* fully analytic computation */
  {
    if (swap == 0) 
      computeKijAnalytic(I0,I1,a,b,c,d);
    else
      computeKijSwappedAnalytic(I0,I1,a,b,c,d);
  }
  else            /* fully analytic or semianalytic */
  {
    if(distanceSegmentToSegment(a,b,c,d)*eta >= sqrt(hi))
    {
      if (swap == 0)
        computeKijSemianalytic(I0,I1,a,b,c,d);
      else
        computeKijSwappedSemianalytic(I0,I1,a,b,c,d);
    }
    else
    {
      if (swap == 0)
        computeKijAnalytic(I0,I1,a,b,c,d);
      else
        computeKijSwappedAnalytic(I0,I1,a,b,c,d);
    }
  }
}


//-----------------------------------------------------------------------------
void computeKijAnalytic(double* I0, double* I1, 
			const Eigen::Vector2d& a, const Eigen::Vector2d& b,
			const Eigen::Vector2d& c, const Eigen::Vector2d& d)
{
  
  double hi = (b-a).squaredNorm(); /* hi = norm(b-a)^2 */
  double hj = (d-c).squaredNorm(); /* hj = norm(d-c)^2 */

  Eigen::Vector2d n = unitNormal(c,d); /* normal vector */
  
  Eigen::Vector2d u = a-b; 
  Eigen::Vector2d v = d-c;
  Eigen::Vector2d w = w=c+d-a-b;

  double dot_u_n = u.dot(n);
  double dot_w_n = w.dot(n);
  double dot_wpu_n = (w+u).dot(n);
  double dot_wmu_n = (w-u).dot(n);

  double det = CrossProd2d(u,v);
  
  double lambda=0.0, mu=0.0;
  if (fabs(det) <= EPS*sqrt(hi*hj))  /* u,v linearly dependent */
  {
    if (fabs(u[0]) > fabs(u[1]))
      mu = v[0]/u[0]; 
    else
      mu = v[1]/u[1];

    *I0 = dot_w_n*( dlp(0,u,w+v)+dlp(0,u,w-v)+ mu*(dlp(1,v,w-u)-dlp(1,v,w+u)) );
    *I1 = dot_w_n*( dlp(0,u,w+v)-dlp(0,u,w-v)+ mu*(dlp(2,v,w-u)-dlp(2,v,w+u)) )*0.5;
  }
  else                               /* u,v linearly independent */
  {
    if (a == d)
    {
      *I0 = 2*( dot_wpu_n*dlp(0,v,w+u)+dot_u_n*dlp(1,u,w-v)+dot_w_n*dlp(0,u,w-v) );
      *I1 =     dot_wpu_n*dlp(1,v,w+u)-dot_u_n*dlp(1,u,w-v)-dot_w_n*dlp(0,u,w-v)
	    + 0.5*(*I0);
    }
    else if (b == c)
    {
      *I0 = 2*( dot_wmu_n*dlp(0,v,w-u)+dot_u_n*dlp(1,u,w+v)+dot_w_n*dlp(0,u,w+v) );
      *I1 =     dot_wmu_n*dlp(1,v,w-u)+dot_u_n*dlp(1,u,w+v)+dot_w_n*dlp(0,u,w+v)
	    - 0.5*(*I0);
    }
    else
    {
      mu     = CrossProd2d(w,v)/det;
      lambda = CrossProd2d(u,w)/det;
     
      *I0 = (mu+1)*dot_wpu_n*dlp(0,v,w+u) - (mu-1)*dot_wmu_n*dlp(0,v,w-u)
            + (lambda+1)*( dot_u_n*dlp(1,u,w+v) + dot_w_n*dlp(0,u,w+v) )
            - (lambda-1)*( dot_u_n*dlp(1,u,w-v) + dot_w_n*dlp(0,u,w-v) );
      *I1 = 0.5*( (mu+1)*dot_wpu_n*dlp(1,v,w+u) - (mu-1)*dot_wmu_n*dlp(1,v,w-u)
                 + (lambda+1)*( dot_u_n*dlp(1,u,w+v) + dot_w_n*dlp(0,u,w+v) )
                 + (lambda-1)*( dot_u_n*dlp(1,u,w-v) + dot_w_n*dlp(0,u,w-v) ) 
                 - lambda*(*I0) ); 
    }
  }
  *I0 *= -0.125*sqrt(hi*hj)/M_PI;
  *I1 *= -0.125*sqrt(hi*hj)/M_PI;
}


//-----------------------------------------------------------------------------
void computeKijSwappedAnalytic(double* I0, double* I1, 
                               const Eigen::Vector2d& a, const Eigen::Vector2d& b,
                               const Eigen::Vector2d& c, const Eigen::Vector2d& d)
{

  double hi = (b-a).squaredNorm(); /* hi = norm(b-a)^2 */
  double hj = (d-c).squaredNorm(); /* hj = norm(d-c)^2 */

  Eigen::Vector2d n; /* normal vector */
  n << (d[1]-c[1])/sqrt(hj) , -(d[0]-c[0])/sqrt(hj);

  Eigen::Vector2d u = d-c;
  Eigen::Vector2d v = a-b;
  Eigen::Vector2d w = w=c+d-a-b;

  double dot_v_n = v.dot(n);       /* dot_v_n=<v,n> */
  double dot_w_n = w.dot(n);
  double dot_wpv_n = (w+u).dot(n); /* dot_wpu_n=<w+u,n> */
  double dot_wmv_n = (w-v).dot(n);

  double det = CrossProd2d(u,v);
  
  double lambda=0.0, mu=0.0;
  if(fabs(det)<=EPS*sqrt(hi*hj))   /* u,v linearly dependent */
  {
    if(fabs(u[0])>fabs(u[1]))
      mu = v[0]/u[0];
    else
      mu = v[1]/u[1];

    *I0 = dot_w_n*( dlp(0,u,w+v)+dlp(0,u,w-v) + mu*( dlp(1,v,w-u)-dlp(1,v,w+u)) );
    *I1 = dot_w_n*( dlp(1,u,w+v)+dlp(1,u,w-v) + mu*(-dlp(1,v,w-u)-dlp(1,v,w+u)
                    +0.5*( dlp(0,u,w+v)-dlp(0,u,w-v)
                           + mu*(dlp(2,v,w-u)-dlp(2,v,w+u)) ) ));
  }
  else                             /* u,v linearly independent */
  {
    if (a == d)
    {
      *I0 = 2*( dot_v_n*dlp(1,v,w-u) + dot_w_n*dlp(0,v,w-u)
                + dot_wpv_n*dlp(0,u,w+v) );
      *I1 =     dot_wpv_n*dlp(1,u,w+v)-dot_v_n*dlp(1,v,w-u)
              - dot_w_n*dlp(0,v,w-u) + 0.5*(*I0);
    }
    else if (b == c)
    {
      *I0 = 2*( dot_v_n*dlp(1,v,w+u) + dot_w_n*dlp(0,v,w+u)
                +dot_wmv_n*dlp(0,u,w-v) );
      *I1 =     dot_v_n*dlp(1,v,w+u) + dot_w_n*dlp(0,v,w+u)
                + dot_wmv_n*dlp(1,u,w-v) - 0.5*(*I0);
    }
    else
    {
      mu     = CrossProd2d(w,v)/det;
      lambda = CrossProd2d(u,w)/det;
     
      *I0 =   (mu+1)*(dot_v_n*dlp(1,v,w+u) + dot_w_n*dlp(0,v,w+u))
            - (mu-1)*(dot_v_n*dlp(1,v,w-u) + dot_w_n*dlp(0,v,w-u))
            + (lambda+1)*dot_wpv_n*dlp(0,u,w+v)
            - (lambda-1)*dot_wmv_n*dlp(0,u,w-v);
      *I1 = 0.5*( (mu+1)*(dot_v_n*dlp(1,v,w+u) + dot_w_n*dlp(0,v,w+u))
                 +(mu-1)*(dot_v_n*dlp(1,v,w-u) + dot_w_n*dlp(0,v,w-u))
                 +(lambda+1)*dot_wpv_n*dlp(1,u,w+v)
                 -(lambda-1)*dot_wmv_n*dlp(1,u,w-v)-mu*(*I0));
    }
  }
  *I0 *= -0.125*sqrt(hi*hj)/M_PI;
  *I1 *= -0.125*sqrt(hi*hj)/M_PI;
}


//-----------------------------------------------------------------------------
void computeKijSemianalytic(double* I0, double* I1, 
                            const Eigen::Vector2d& a, const Eigen::Vector2d& b,
                            const Eigen::Vector2d& c, const Eigen::Vector2d& d)
{
  double hi = (b-a).squaredNorm(); /* hi = norm(b-a)^2 */
  double hj = (d-c).squaredNorm(); /* hj = norm(d-c)^2 */
  
  /* 16-point Gaussian quadrature on [-1,1] */
  const int order = 16;
  const double* gauss_weight = getGaussWeights(order);
  const double* gauss_point  = getGaussPoints(order);

  Eigen::Vector2d n = unitNormal(c,d);

  Eigen::Vector2d u = a-b;
  Eigen::Vector2d v = d-c;
  Eigen::Vector2d w = w=c+d-a-b;

  double I0tmp = 0.0;
  double I1tmp = 0.0;
  for (int i=0;i<order;++i)
  {
    Eigen::Vector2d z = gauss_point[i]*u + w;
   
    double dot_z_n = z.dot(n);
   
    I0tmp += gauss_weight[i]*dot_z_n*dlp(0,v,z);
    I1tmp += gauss_weight[i]*dot_z_n*dlp(1,v,z);
  }

  *I0 = -0.125*sqrt(hi*hj)*I0tmp/M_PI;
  *I1 = -0.125*sqrt(hi*hj)*I1tmp/M_PI;
}


//-----------------------------------------------------------------------------
void computeKijSwappedSemianalytic(double* I0, double* I1, 
                                   const Eigen::Vector2d& a,
                                   const Eigen::Vector2d& b,
                                   const Eigen::Vector2d& c,
                                   const Eigen::Vector2d& d)
{
  double hi = (b-a).squaredNorm(); /* hi = norm(b-a)^2 */
  double hj = (d-c).squaredNorm(); /* hj = norm(d-c)^2 */
  
  /* 16-point Gaussian quadrature on [-1,1] */
  const int order = 16;
  const double* gauss_weight = getGaussWeights(order);
  const double* gauss_point  = getGaussPoints(order);

  Eigen::Vector2d n = unitNormal(c,d); /* normal vector */

  Eigen::Vector2d u = d-c;
  Eigen::Vector2d v = a-b;
  Eigen::Vector2d w = w=c+d-a-b;

  double dot_v_n = v.dot(n);
  double dot_w_n = w.dot(n);
 
  double I0tmp = 0.0;
  double I1tmp = 0.0;
  for (int i=0;i<order;++i)
  {
    Eigen::Vector2d z = gauss_point[i]*u + w;
   
    I0tmp += gauss_weight[i]*( dot_v_n*dlp(1,v,z) + dot_w_n*dlp(0,v,z) );
    I1tmp += gauss_weight[i]*gauss_point[i]*( dot_v_n*dlp(1,v,z) 
                                             + dot_w_n*dlp(0,v,z) );
  }

  *I0 = -0.125*sqrt(hi*hj)*I0tmp/M_PI;
  *I1 = -0.125*sqrt(hi*hj)*I1tmp/M_PI;
}

