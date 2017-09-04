///////////////////////////////////////////////////////////////////////////////
/// \file newtonPotential.cpp
/// \brief This file provides functions to compute the newton potential.
///
///  This file is part of the HILBERT program package for the numerical
///  solution of the Laplace equation with mixed boundary conditions by use of
///  BEM in 2D.
///
///  C++ adaptation for ANCSE17 of HILBERT V3.1 TUWien 2009-2013
///////////////////////////////////////////////////////////////////////////////
#include <cstdlib>

#include "constants.hpp"
#include "newtonPotential.hpp"

void computeN(double* N, const double *coordinates, const double *elements,
    const double* vertices, const double* triangles, double eta) {

  int j=0, k=0;
  int aidx=0, bidx=0, n1idx=0, n2idx=0, n3idx=0;
  double a[2], b[2], n1[2], n2[2], n3[2];

  for (j=0;j<nT;++j) { /*  running over triangles */
    n1idx = (int) triangles[j]-1;
    n2idx = (int) triangles[j+nT]-1;
    n3idx = (int) triangles[j+2*nT]-1;

    n1[0] = vertices[n1idx];
    n1[1] = vertices[n1idx+nV];
    n2[0] = vertices[n2idx];
    n2[1] = vertices[n2idx+nV];
    n3[0] = vertices[n3idx];
    n3[1] = vertices[n3idx+nV];

    for (k=0;k<nE;++k) { /* running over boundary elements */
      aidx = (int) elements[k]-1;
      bidx = (int) elements[k+nE]-1;

      a[0] = coordinates[aidx];
      a[1] = coordinates[aidx+nC];
      b[0] = coordinates[bidx];
      b[1] = coordinates[bidx+nC];

      N[k+j*nE] = computeNkj(a,b,n1,n2,n3,eta);
    }
  }
}

double computeNkj(const Eigen::Matrix2d &a, const Eigen::Matrix2d &b,
                  const Eigen::MatrixXd &nodes, double eta)
{
  double diam = 0.;
  double tmp = 0.;
  double x0 = 0., x1 = 0.;

  double distET = distanceSegmentToSegment(a[0], a[1], b[0], b[1],
                                          nodes[0], nodes[1], n2[0], n2[1]);
  tmp = distanceSegmentToSegment(a[0],a[1],b[0],b[1],
                                  nodes[0],nodes[1],n3[0],n3[1]);
  if (tmp < distET)
    distET = tmp;
  
  tmp = distanceSegmentToSegment(a[0],a[1],b[0],b[1],
                                  n2[0],n2[1],n3[0],n3[1]);
  if (tmp < distET)
    distET = tmp;
  
  x0 = a[0] - b[0]; x1 = a[1] - b[1];
  diam = x0*x0+x1*x1;
 
  if (distET*eta >= diam) 	/* Semi-analytic */
     return computeNkjSemiAnalyticSegment(a,b,nodes,n2,n3);
  else 				/* analytic */
      return computeNkjAnalytic(a,b,nodes,n2,n3);
}

double computeNkjSemiAnalyticSegment(const Eigen::Matrix2d& a,
                                     const Eigen::Matrix2d& b,
                                     const Eigen::MatrixXd& nodes) {

  int k;
  const double *gauss_points, *gauss_weights;
  double tmp[2], val;

  gauss_points = getGaussPoints(GAUSS_ORDER);
  gauss_weights = getGaussWeights(GAUSS_ORDER);

  val=0.;
  for (k=0;k<GAUSS_ORDER;++k){
    tmp[0] = 0.5*(a[0]+b[0] + gauss_points[k]*(b[0]-a[0]));
    tmp[1] = 0.5*(a[1]+b[1] + gauss_points[k]*(b[1]-a[1]));

    val = val + gauss_weights[k]*newtonPotential(n1,n2,n3,tmp);
  }
  return -val*sqrt((a[0]-b[0])*(a[0]-b[0])+(a[1]-b[1])*(a[1]-b[1]))*
            0.25 / M_PI;
}

double computeNkjSemiAnalyticTriangle(double a[2], double b[2],
                                      double n1[2], double n2[2],
                                      double n3[2])
{

  double detFT = 0.;
  const int nr_gauss_points = 7;
  const double* gauss_point_x;
  const double* gauss_point_y;
  const double* gauss_weights;
  double aux1[2], aux2[2];
  double tmp;
  int i;

  gauss_point_x = getGaussPointsT(nr_gauss_points,0);
  gauss_point_y = getGaussPointsT(nr_gauss_points,1);
  gauss_weights = getGaussWeightsT(nr_gauss_points);

  detFT = fabs((n2[0]-n1[0])*(n3[1]-n1[1])-(n3[0]-n1[0])*(n2[1]-n1[1]));
  tmp=0.;

  for (i=0;i<nr_gauss_points;i++)
  {
    aux1[0] = 0.5*(b[0]-a[0]);
    aux1[1] = 0.5*(b[1]-a[1]);
    aux2[0] = 0.5*(b[0]+a[0])-(n2[0]-n1[0])*gauss_point_x[i] - (n3[0]-n1[0])
                                                  *gauss_point_y[i] - n1[0];
    aux2[1] = 0.5*(b[1]+a[1])-(n2[1]-n1[1])*gauss_point_x[i] - (n3[1]-n1[1])
                                                  *gauss_point_y[i] - n1[1];
    tmp = tmp + gauss_weights[i]*slp(0,aux1,aux2);
  }

  return -0.125*sqrt((a[0]-b[0])*(a[0]-b[0])+(a[1]-b[1])*(a[1]-b[1]))*detFT*
            tmp/M_PI;
}


double computeNkjAnalytic(const Eigen::Matrix2d &a, const Eigen::Matrix2d &b,
                          const Eigen::MatrixXd &nodes) {

  double detFK = 0., volK = 0., lengthT = 0.;
  double normasquared;
  double u[2], n1minusn2over2a[2], n1minusn3over2[2], mj[2],
          aux[2], aux1[2], aux2[2], bjminusajover2[2];
  double n2minusn3over2[2];
  double coeff1=0., coeff2=0., coeff3=0., coeff4=0., coeff5=0., coeff6=0.;
  double J1=0., J2=0., J3=0., J4=0., J5=0., J6=0.;
  double atanInt=0., ret=0.;

  detFK = fabs((n2[0]-nodes[0])*(n3[1]-nodes[1])-(n3[0]-nodes[0])*(n2[1]-nodes[1]));
  volK = 0.5*detFK;
  lengthT = sqrt((a[0]-b[0])*(a[0]-b[0])+(a[1]-b[1])*(a[1]-b[1]));

  u[0] = nodes[0]-n2[0]; u[1] = nodes[1]-n2[1];
  normasquared = dot(u,u);
  n1minusn2over2a[0] = 0.5*u[0]/normasquared;
  n1minusn2over2a[1] = 0.5*u[1]/normasquared;
  n1minusn3over2[0] = 0.5*(nodes[0]-n3[0]);
  n1minusn3over2[1] = 0.5*(nodes[1]-n3[1]);
  n2minusn3over2[0] = 0.5*(n2[0]-n3[0]);
  n2minusn3over2[1] = 0.5*(n2[1]-n3[1]);

  /* aux holds (n1-n3)/2 -n1 + mj */
  mj[0] = 0.5*(a[0]+b[0]); mj[1] = 0.5*(a[1]+b[1]);
  aux[0] = n1minusn3over2[0] - nodes[0] + mj[0];
  aux[1] = n1minusn3over2[1] - nodes[1] + mj[1];

  /* aux1 holds mj - (n2+n3)/2 */
  aux1[0] = mj[0] - 0.5*(n2[0]+n3[0]); aux1[1] = mj[1] - 0.5*(n2[1]+n3[1]);

  /* aux2 holds mj - (n1+n3)/2 */
  aux2[0] = mj[0] - 0.5*(nodes[0]+n3[0]); aux2[1] = mj[1] - 0.5*(nodes[1]+n3[1]);

  bjminusajover2[0] = 0.5*(b[0]-a[0]); bjminusajover2[1] = 0.5*(b[1]-a[1]);

  coeff1 = 0.25 + dot(n1minusn2over2a,aux);
  coeff2 = -0.25 + dot(n1minusn2over2a,n1minusn3over2);
  coeff3 = dot(n1minusn2over2a,bjminusajover2);
  coeff4 = dot(n1minusn2over2a,aux);
  coeff5 = dot(n1minusn2over2a,n1minusn3over2);
  coeff6 = coeff3;

  J1 = doubleSlp(0,0,bjminusajover2,n2minusn3over2,aux1);
  J2 = doubleSlp(0,1,bjminusajover2,n2minusn3over2,aux1);
  J3 = doubleSlp(1,0,bjminusajover2,n2minusn3over2,aux1);
  J4 = doubleSlp(0,0,bjminusajover2,n1minusn3over2,aux2);
  J5 = doubleSlp(0,1,bjminusajover2,n1minusn3over2,aux2);
  J6 = doubleSlp(1,0,bjminusajover2,n1minusn3over2,aux2);

  atanInt = integrateAtanInt(a,b,nodes,n2,n3);
  ret = lengthT*(coeff1*J1 + coeff2*J2 + coeff3*J3 - coeff4*J4
                 - coeff5*J5 - coeff6*J6 - 1)
        + 0.5*atanInt/(normasquared);

  return -0.5*volK*ret/M_PI;
}

double newtonPotential(const Eigen::MatrixXd &nodes, const Eigen::Matrix2d &x)
{
  double u[2], v[2], w[2];
  double a = 0;
  double atanint=0.;
  double result = 0.;
  double integralI,integralII,integralIII;
  
  double v_minus_u_div_2[2];
  double w_plus_u_plus_v_minus_u_div_2[2];
  double v_div_2[2];
  double w_plus_v_div_2[2];
  
  double* slpA = NULL;
  double* slpB = NULL;
  
  double volT = 0.5*fabs(
        (n2[0]-nodes[0])*(n3[1]-nodes[1])-(n3[0]-nodes[0])*(n2[1]-nodes[1]));

  u[0] = nodes[0]-n2[0];
  u[1] = nodes[1]-n2[1];
  
  v[0] = nodes[0]-n3[0];
  v[1] = nodes[1]-n3[1];
  
  w[0] = x[0]-nodes[0];
  w[1] = x[1]-nodes[1];
  
  a = u[0]*u[0]+u[1]*u[1];
  
  v_minus_u_div_2[0] = (v[0]-u[0])/2.;
  v_minus_u_div_2[1] = (v[1]-u[1])/2.;
  
  w_plus_u_plus_v_minus_u_div_2[0] = w[0]+u[0]+v_minus_u_div_2[0];
  w_plus_u_plus_v_minus_u_div_2[1] = w[1]+u[1]+v_minus_u_div_2[1];
  
  v_div_2[0] = v[0]/2.;
  v_div_2[1] = v[1]/2.;
  
  w_plus_v_div_2[0] = w[0] + v_div_2[0];
  w_plus_v_div_2[1] = w[1] + v_div_2[1];
  
  slpA = slpIterative(1,v_minus_u_div_2,w_plus_u_plus_v_minus_u_div_2);
  slpB = slpIterative(1,v_div_2,w_plus_v_div_2);
  
  integralI = 0.25*(slpA[0]-slpA[1]);
  integralII = (u[0]*w[0]+u[1]*w[1])*slpA[0]+
               0.5*(u[0]*v[0]+u[1]*v[1])*slpA[0]+
               0.5*(u[0]*v[0]+u[1]*v[1])*slpA[1];
  integralIII = (u[0]*w_plus_v_div_2[0]+u[1]*w_plus_v_div_2[1])*slpB[0]+
                +0.5*(u[0]*v[0]+u[1]*v[1])*slpB[1];
  atanint = evalAtanInt(nodes,n2,n3,x);

  result = volT*integralI + volT/(2.*a)*integralII -
           volT/(2.*a)*integralIII - volT + volT/(2*a)*atanint;
  
  free(slpA); free(slpB);

  return result;
}

double evalAtanInt(const Eigen::MatrixXd& nodes, const Eigen::Matrix2d& x) {

  double u[2], v[2], w[2];
  const double* gauss_point;
  const double* gauss_wht;
  double val=0., xi;
  double valInner=0., uu, vv, u_dot_v;
  double a,b,c;
  int j;

  gauss_point = getGaussPoints(GAUSS_ORDER);
  gauss_wht   = getGaussWeights(GAUSS_ORDER);

  /* compute u and v */
  u[0] = n1[0] - n2[0];
  u[1] = n1[1] - n2[1];
  v[0] = n1[0] - n3[0];
  v[1] = n1[1] - n3[1];
  uu = u[0]*u[0] + u[1]*u[1];
  vv = v[0]*v[0] + v[1]*v[1];
  u_dot_v = u[0]*v[0] + u[1]*v[1];

  /* compute w */
  w[0] = x[0] - n1[0];
  w[1] = x[1] - n1[1];

  val=0.;

  for ( j= 0; j<GAUSS_ORDER; ++j ){
    /* transformation of quadrature points from [-1,1] to [0,1] */
    xi = (gauss_point[j]+1)*0.5;

    a = uu;
    b = 2*(u[0]*w[0]+u[1]*w[1]+xi*(u[0]*v[0]+u[1]*v[1]));
    c = w[0]*w[0]+w[1]*w[1]+xi*xi*vv+2*xi*(w[0]*v[0]+w[1]*v[1]);

    valInner = innerAtanInt(a,b,c,xi);
    val += gauss_wht[j]*valInner;
  }

  return val;
}

double innerAtanInt(double a, double b, double c, double xi)
{
  double delta = 4*a*c-b*b;
  double deltaRoot = sqrt(delta);
  double s = (2*a*(1-xi)+b)/deltaRoot;
  double t = -b/deltaRoot;
  double st = s*t;
  double valInner = 0.;

  if(delta>EPS*1*a*c) {
    if ( st<1-EPS )
      valInner = deltaRoot*atan( (s+t)/(1-st));
    else {
      if ( st>1+EPS) {
        if (s < 0)
	  valInner = deltaRoot*(atan( (s+t)/(1-st)) - M_PI);
        else
	  valInner = deltaRoot*(atan( (s+t)/(1-st)) + M_PI);
      }
      else {
        valInner = deltaRoot*(atan(s) + atan(t));
      }
    }
  }
  return valInner;
}


double integrateAtanInt(const Eigen::Matrix2d &a, const Eigen::Matrix2d &b,
                        const Eigen::MatrixXd &nodes)
{
  double u[2], v[2];
  const double* gauss_point;
  const double* gauss_wht;
  int k;
  double sx[2], val=0.;
  double valOuter=0., uu, vv, uxv;

  gauss_point = getGaussPoints(GAUSS_ORDER);
  gauss_wht   = getGaussWeights(GAUSS_ORDER);

  /* compute u and v */
  u[0] = nodes[0] - n2[0];
  u[1] = nodes[1] - n2[1];
  v[0] = nodes[0] - n3[0];
  v[1] = nodes[1] - n3[1];
  uu = u[0]*u[0] + u[1]*u[1];
  vv = v[0]*v[0] + v[1]*v[1];
  uxv = u[0]*v[0] + u[1]*v[1];

  val=0.;

  for ( k=0; k<GAUSS_ORDER; ++k )
  {
    /* transformation of quadrature points from [-1,1] to [a,b] */
    sx[0] = ((1-gauss_point[k])*a[0]+(1+gauss_point[k])*b[0])*0.5;
    sx[1] = ((1-gauss_point[k])*a[1]+(1+gauss_point[k])*b[1])*0.5;

    valOuter = evalAtanInt(nodes,n2,n3,sx);

    val += gauss_wht[k]*valOuter;
  }
  
  val = val*0.5*sqrt((a[0]-b[0])*(a[0]-b[0]) + (a[1]-b[1])*(a[1]-b[1]));
  return val;
}

