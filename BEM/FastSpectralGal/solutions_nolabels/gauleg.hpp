//// 
//// Copyright (C) 2017 SAM (D-MATH) @ ETH Zurich
//// Author(s): curzuato < > 
//// Contributors:  dcasati 
//// This file is part of the AdvNumCSE repository.
////
#ifndef GAULEGHPP
#define GAULEGHPP

#include <iostream>
#include <cmath>
#include <Eigen/Dense>


/* @brief Compute Gaussian quadrature nodes and weights for n nodes over interval [a,b]
 * \param[in] a,b Interval [a,b] endpoints
 * \param[in] n Number of quadrature points
 * \param[out] xq,wq Pair of Gauss-Legendre quadrature points and weights
 */
inline std::pair<Eigen::RowVectorXd,Eigen::RowVectorXd>
gauleg(double a,double b,int n, double eps=1.e-13) {

  const double PI = 4.*std::atan(1.); // PI
  int i,j,m;
  double xmid,xlen,p1,p2,p3,dp1,z,z1,wqi;
  Eigen::RowVectorXd xq(n),wq(n);

  m = (n + 1)/2;
  xmid = 0.5*(a + b);
  xlen = 0.5*(b - a);

   // get roots
   for (i = 0; i < m; i++){

     // i-th root guess
     z = std::cos(PI*(i + 1 - 0.25)/(n + 0.5));

     // get i-th root
     do {
       p1 = 1.;
       p2 = 0.;
       for (j = 1; j <= n; j++) {
         p3 = p2;
         p2 = p1;
         p1 = ((2.*j - 1.)*z*p2 - (j - 1.)*p3)/j;
       }
       dp1 = n*(z*p1 - p2)/(z*z - 1.0);
       z1 = z;
       z = z1 - p1/dp1;
     } while ( std::abs(z - z1) > eps);

     // set nodes
     xq(    i) = xmid - xlen*z;
     xq(n-1-i) = xmid + xlen*z;

     // set weights
     wqi = 2.*xlen/((1. - z*z)*dp1*dp1);
     wq(    i) = wqi;
     wq(n-1-i) = wqi;
   }

  // return
  return std::make_pair(xq,wq);

}

#endif
