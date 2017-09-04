///////////////////////////////////////////////////////////////////////////////
/// \file singleLayerPotential.cpp
/// \brief This file provides functions to compute single layer potential (slp)
///  type integrals.
///
///  This file contains only the implementation. For detailed documentation see
///  singleLayerPotential.hpp
///
///  This file is part of the HILBERT program package for the numerical
///  solution of the Laplace equation with mixed boundary conditions by use of
///  BEM in 2D.
///
///  C++ adaptation for ANCSE17 of HILBERT V3.1 TUWien 2009-2013
///////////////////////////////////////////////////////////////////////////////

#include <cstdlib>

#include "singleLayerPotential.hpp"
#include "constants.hpp"


//------------------------------------------------------------------------------
double slp(int k, const Eigen::Vector2d& u, const Eigen::Vector2d& v)
{
  Eigen::VectorXd tmp = slpIterative(k, u, v);
  return tmp[k];
}

//------------------------------------------------------------------------------
Eigen::VectorXd slpIterative(int k, const Eigen::Vector2d& u,
			     const Eigen::Vector2d& v)
{
  double a = u.squaredNorm();  /* a = <u,u> */
  double b = 2 * u.dot(v);     /* b = 2 <u,v> */
  double c = v.squaredNorm();  /* c = <v,v> */
  double D = 0.;
  Eigen::VectorXd val(k+1);

  /* Ensure that discriminant is either positive or zero */
  double tmp = 4*a*c - b*b;
  assert(fabs(u[0]) > EPS || fabs(u[1]) > EPS
          || fabs(v[0]) > EPS || fabs(v[1]) > EPS);
  assert(tmp >= -fabs(EPS*4*a*c)); /* By theory there holds tmp >= 0. */

  if (tmp > EPS*4*a*c)
    D = sqrt(tmp);
  else
    D = 0.;
  
  /* The case k=0 */
  if (fabs(u[0]) < EPS && fabs(u[1]) < EPS) {
      val[0] = 2*log(c);
  }
  else if (D == 0.) {
    tmp = b + 2*a;
    if (fabs(tmp) > EPS*a)
      val[0] = tmp * log( 0.25*tmp*tmp /a );
    else
      val[0] = 0;
    tmp = b - 2*a;
    if (fabs(tmp) > EPS*a)
      val[0] -= tmp * log( 0.25*tmp*tmp /a );
    val[0] = 0.5*val[0] /a - 4;
  }
  else { /* case D > 0 */
    tmp = c - a;
    if (fabs(tmp) < EPS*c)
      val[0] = 0.5*M_PI;
    else if (a < c)
      val[0] = atan( D /tmp );
    else
      val[0] = atan( D /tmp ) + M_PI;

    val[0] = ( 0.5*( (b+2*a) * log(a+b+c) - (b-2*a) * log(a-b+c) )
                + D*val[0]) / a - 4;
  }
  if (k == 0)
    return val;

  /* The case k=1 */
  if (k>=1) {
    if (fabs(u[0]) < EPS && fabs(u[1]) < EPS) {
      val[1] = 0.;
    }
    else {
      /* val holds \int_{-1}^{+1} \log |a*s^2+b*s+c|^2 ds. */
      val[1] = -b*(2+val[0]);

      tmp = a+b+c;
      if (fabs(tmp) > EPS*a)
        val[1] += tmp * log(tmp);
   
      tmp = a-b+c;
      if (fabs(tmp) > EPS*a)
        val[1] -= tmp * log(tmp);

      val[1] /= (2*a);
    }
  }
  if (k == 1)
    return val;

  /* The case k>=2 */
  for (int i=2; i <= k; ++i) {
    if (fabs(u[0]) < EPS && fabs(u[1]) < EPS) {
      if (i%2 == 0)
        val[i] = 2.*log(c)/(i+1);
      else
        val[i] = 0.;
    }
    else {
      tmp = a+b+c;
      if (fabs(tmp) > a*EPS)
        val[i] = tmp*log(tmp);
      else
        val[i] = 0.;

      tmp = a-b+c;
      if (i % 2 == 0) {
        if (fabs(tmp) > a*EPS)
          val[i] += tmp*log(tmp);
        val[i] -= 4*a/(i+1);
      }
      else {
        if (fabs(tmp) > a*EPS)
          val[i] -= tmp*log(tmp);
        val[i] -= 2*b/i;
      }

      val[i] -= i*b*val[i-1]+(i-1)*c*val[i-2];
      val[i] /= ((i+1)*a);
    }
  }

  return val;
}

//------------------------------------------------------------------------------
double doubleSlp(int k, int l, const Eigen::Vector2d& u,
		 const Eigen::Vector2d& v, const Eigen::Vector2d& w)
{

  double output = 0.;
  Eigen::VectorXd memTableUWpv, memTableUWmv, memTableVWpu, memTableVWmu;
  double normUSq = u.squaredNorm();
  double normVSq = v.squaredNorm();
  double normWSq = w.squaredNorm();
  double detUV = u[0]*v[1] - u[1]*v[0];

  if (normUSq < EPS && normVSq < EPS)
  {
    if (k%2 == 0)
    {
      if (l%2 == 0)
      {
        if (normWSq < EPS)
          return 0;
        else
          return 2./((k+1)*(l+1))*log(w[0]*w[0]+w[1]*w[1]);
      }
    }
    return 0.;
  }
  else if (normUSq < EPS)
  {
    if (k%2 == 0)
      return 1./(k+1) * slp(l, v, w);
    else
      return 0.;
  }
  else if (normVSq < EPS)
  {
    if (l%2 == 0)
      return 1./(l+1) * slp(k, u, w);
    else
      return 0.;
  }

  if (fabs(detUV) < EPS*sqrt(normUSq*normVSq)) { /* u,v parallel */
    double mu = 0.;

    memTableUWpv = slpIterative(k, u, w+v);
    memTableUWmv = slpIterative(k, u, w-v);
    memTableVWpu = slpIterative(k+l+1, v, w+u);
    memTableVWmu = slpIterative(k+l+1, v, w-u);

    if (fabs(u[0]) < fabs(u[1]))
      mu = v[1] / u[1];
    else
      mu = v[0] / u[0];

    output = memTableUWpv[0] - mu*memTableVWpu[k+l+1] +
              mu*memTableVWmu[k+l+1];
    if ((k+l) % 2 == 0)
      output += memTableUWmv[0];
    else
      output -= memTableUWmv[0];
    output /= (2*(k+l+1));

    for (int i = 1; i <= k; ++i) {
      output *= 2*i*mu;
      output += memTableUWpv[i] - mu*memTableVWpu[l+k+1-i]; 

      if ((k+l-i) % 2 == 0)
        output += memTableUWmv[i];
      else
        output -= memTableUWmv[i];

      if (i % 2 == 0)
        output += mu*memTableVWmu[l+k+1-i];
      else
        output -= mu*memTableVWmu[l+k+1-i];

      output /= (2*(l+k+1-i));
    }
  }
  else {
    double mu1 = 0., mu2 = 0.;
    Eigen::VectorXd tmp(l+1);

    if (fabs(w[0]+v[0]-u[0]) < fabs(w[0])*EPS
          && fabs(w[1]+v[1]-u[1]) < fabs(w[1])*EPS) {
      memTableVWmu.resize(l+1);
      for (int i = 0; i <= l; ++i)
        memTableVWmu[i] = 0.;

      memTableUWpv.resize(k+1);
      for (int i = 0; i <= k; ++i)
        memTableUWpv[i] = 0.;
    }
    else {
      memTableVWmu = slpIterative(l, v, w-u);
      memTableUWpv = slpIterative(k, u, w+v);
    }

    if (fabs(w[0]+u[0]-v[0]) < fabs(w[0])*EPS
          && fabs(w[1]+u[1]-v[1]) < fabs(w[1])*EPS) {
      memTableVWpu.resize(l+1);
      for (int i = 0; i <= l; ++i)
        memTableVWpu[i] = 0.;

      memTableUWmv.resize(k+1);
      for (int i = 0; i <= k; ++i)
        memTableUWmv[i] = 0.;
    }
    else {
      memTableVWpu = slpIterative(l, v, w+u);
      memTableUWmv = slpIterative(k, u, w-v);
    }

    mu1 = ( v[1]*w[0] - v[0]*w[1]) / detUV;
    mu2 = (-u[1]*w[0] + u[0]*w[1]) / detUV;

    tmp[0] = -2 + ((mu1+1)*memTableVWpu[0] - (mu1-1)*memTableVWmu[0]
              + (mu2+1)*memTableUWpv[0] - (mu2-1)*memTableUWmv[0]) * 0.25;

    for (int i = 1; i <= l; ++i) {
      tmp[i] = 0.5*((mu1+1)*memTableVWpu[i] - (mu1-1)*memTableVWmu[i]
            + (mu2+1)*memTableUWpv[0]) - i*mu2*tmp[i-1];
      if (i%2 == 0) {
        tmp[i] -= 4./(i+1);
        tmp[i] -= 0.5 * (mu2-1)*memTableUWmv[0];
      }
      else
        tmp[i] += 0.5 * (mu2-1)*memTableUWmv[0];

      tmp[i] /= (i+2);
    }

    for (int i = 1; i <= k; ++i) {
      tmp[0] = 0.5*((mu1+1)*memTableVWpu[0] + (mu2+1)*memTableUWpv[i]
                  - (mu2-1)*memTableUWmv[i]) - i*mu1*tmp[0];
      if (i%2 == 0) {
        tmp[0] -= 4./(i+1);
        tmp[0] -= 0.5*(mu1-1)*memTableVWmu[0];
      }
      else {
        tmp[0] += 0.5*(mu1-1)*memTableVWmu[0];
      }

      tmp[0] /= (i+2);

      for (int j = 1; j <= l; ++j) {
        tmp[j] *= -i*mu1;
        tmp[j] -= j*mu2*tmp[j-1];
        tmp[j] += 0.5*( (mu1+1)*memTableVWpu[j] + (mu2+1)*memTableUWpv[i] );
        if (i%2 == 0) {
          if (j%2 == 0) {
            tmp[j] -= 4./((i+1)*(j+1));
          }
          tmp[j] -= 0.5 * (mu1-1) * memTableVWmu[j];
        }
        else {
          tmp[j] += 0.5 * (mu1-1) * memTableVWmu[j];
        }

        if (j%2 == 0) {
          tmp[j] -= 0.5 * (mu2-1) * memTableUWmv[i];
        }
        else {
          tmp[j] += 0.5 * (mu2-1) * memTableUWmv[i];
        }

        tmp[j] /= (i+j+2);
      }
    }

    output = tmp[l];
  }

  return output;
}

//------------------------------------------------------------------------------
double computeVij(const Eigen::Vector2d& a, const Eigen::Vector2d& b,
		  const Eigen::Vector2d& c, const Eigen::Vector2d& d, double eta)
{
  double hi = (b-a).squaredNorm(); /* hi = norm(b-a)^2 */
  double hj = (d-c).squaredNorm(); /* hj = norm(d-c)^2 */

  return sqrt(hi*hj)*computeWij(a,b,c,d, eta);
}

//------------------------------------------------------------------------------
double computeWij(Eigen::Vector2d a, Eigen::Vector2d b,
		  Eigen::Vector2d c, Eigen::Vector2d d, double eta)
{
  double hi = (b-a).squaredNorm(); /* hi = norm(b-a)^2 */
  double hj = (d-c).squaredNorm(); /* hj = norm(d-c)^2 */
  double tmp = 0.;

  /* For stability reasons, we guarantee   hj <= hi   to ensure that *
   * outer integration is over smaller domain. This is done by       *
   * swapping Ej and Ei if necessary.                                */
  if (hj > hi) {
    std::swap(a,c);
    std::swap(b,d);
    std::swap(hi,hj);
  }

  if ( eta == 0) { /* compute all matrix entries analytically */
    return computeWijAnalytic(a,b,c,d);
  }
  else { /* compute admissible matrix entries semi-analytically */
    if ( distanceSegmentToSegment(a,b,c,d) > eta*sqrt(hj) )
    {
      return computeWijSemianalytic(a,b,c,d);
    }
    else {
      return computeWijAnalytic(a,b,c,d);
    }
  }
}

//------------------------------------------------------------------------------
double computeWijAnalytic(const Eigen::Vector2d& a, const Eigen::Vector2d& b,
			  const Eigen::Vector2d& c, const Eigen::Vector2d& d)
{
  double hi = (b-a).squaredNorm(); /* hi = norm(b-a)^2 */
  double hj = (d-c).squaredNorm(); /* hj = norm(d-c)^2 */
  double val = 0.;
  double lambda, mu;

  Eigen::Vector2d x = (b-a)/2.;
  Eigen::Vector2d y = (c-d)/2.;
  Eigen::Vector2d z = (a+b-c-d)/2.;

  /* There hold different recursion formulae if Ei and Ej */
  /* are parallel (det = 0) or not                        */
  double det = x[0]*y[1] - x[1]*y[0];

  if ( fabs(det) <= EPS*sqrt(hi*hj) ) { /* case that x and y are linearly */
    if ( fabs(x[0]) < fabs(x[1]) )      /* dependent, i.e., Ei and Ej are */
      lambda = y[1] / x[1];             /* parallel. */
    else
      lambda = y[0] / x[0];

    val = 0.5*( lambda * ( slp(1, y, z-x) - slp(1, y, z+x) )
                         + slp(0, x, z+y) + slp(0, x, z-y) );
  }

  else { /* case that x and y are linearly independent */
    lambda = (z[0]*y[1] - z[1]*y[0]) /det;
    mu = (x[0]*z[1] - x[1]*z[0]) /det;

    val = 0.25 * (-8 + (lambda+1)*slp(0, y, z+x) - (lambda-1)*slp(0, y, z-x)
                          + (mu+1)*slp(0, x, z+y) - (mu-1)*slp(0, x, z-y));
  }
  
  return -0.125*val /M_PI; /* = -1/(8*M_PI)*val */
}

//------------------------------------------------------------------------------
double computeWijSemianalytic(const Eigen::Vector2d& a, const Eigen::Vector2d& b,
			      const Eigen::Vector2d& c, const Eigen::Vector2d& d)
{
  double val = 0;
  const double* gauss_point = getGaussPoints(GAUSS_ORDER);
  const double* gauss_wht = getGaussWeights(GAUSS_ORDER);

  Eigen::Vector2d u = 0.5*(a-b);

  for (int k=0; k<GAUSS_ORDER; ++k){
      /* transformation of quadrature nodes from [-1,1] to [a,b] */
      Eigen::Vector2d sx = ((1-gauss_point[k])*c+(1+gauss_point[k])*d)*0.5;      
      Eigen::Vector2d v = sx - 0.5*(a+b);
 
      /* inner product wht*func(sx) */
      val += gauss_wht[k] * slp(0, u, v);
  }
  
  return -0.0625*val / M_PI; /* = - 1/(16*M_PI) * int(log |.|^2) */
}
