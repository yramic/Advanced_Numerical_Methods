//// 
//// Copyright (C) 2017 SAM (D-MATH) @ ETH Zurich
//// Author(s): ascapin < > 
//// Contributors:  ppanchal 
//// This file is part of the AdvNumCSE repository.
////
#include <iostream>
#include <cmath>
#include <Eigen/Dense>
#include "gauleg.hpp"

// Quadrature rule struct
struct QuadRule {
  std::size_t     dim; // dimension of space
  std::size_t     n;   // number of nodes/weights
  Eigen::MatrixXd x;   // quadrature nodes (columns of a matrix with dim rows)
  Eigen::VectorXd w;   // vector of quadrature weights
};


/* @brief Build a Gauss-Legendre quadrature rule on the unit interval
 * \param[in] n Number of quadrature points and weights
 * \param[out] Quadrule Corresponding Gauss-Legendre quadrature rule
 */
QuadRule GauLeg(int n) {

  // declare
  QuadRule QR1D;

    // TODO: Use gauleg to set the quadrature

  // return
  return QR1D;

}


/* @brief Build a quadrature rule on the unit quadrilateral by tensor
 *        product of one-dimensional Gauss-Legendre quadrature rule
 * \param[in] QR1D One-dimendional quadrature rule.
 * \param[out] QuadRule Two-dimensional tensor-product quadrature rule.
 */
QuadRule Quad_QR(QuadRule QR1D) {

  // declare
  QuadRule QR2D;

    // TODO: Use QR1D to set the 2D-quadrature by means of tensor quadrature.
    
  // return
  return QR2D;

};


/* @brief Build a quadrature rule on the unit triangle from a quadrature
 *        rule on the unit triangle
 * \param[in] QR1D One-dimendional quadrature rule.
 * \param[out] QuadRule Two-dimensional quadrature rule constructed with the 
 *             duffy trick.
 */
QuadRule Tria_QR(QuadRule QR1D) {
  // declare
  QuadRule QR2D;

    // TODO: Use Duffy trick to set QR2D from QR1D

  // return
  return QR2D;

};


/* @brief integrate function with given quadrature rule
 * \param[in] QR Quadrature rule
 * \param[in] func Function to be integrated. Dimension of QR and func should match!
 * \param[out] double Obtained value
 */
template <typename Function>
double integrate(QuadRule QR,Function func) {

  double sum = 0.;
    // TODO: Perform numerical integration of func using QR
    
  // return
  return sum;

}

// Purpose: compute factorial (n!)
//
int factorial(int n) {
  return n == 0 ? 1 : n*factorial(n-1);
}


int main() {

  // set number of Gauss-Legendre quadrature nodes
  int N = 2;

  // build quadrature rules
  auto QR1D = GauLeg(N);
  auto Quad_QR2D = Quad_QR(QR1D);
  auto Tria_QR2D = Tria_QR(QR1D);

  // test 1D quadrature on monomials in interval [0,1]
  std::cout << std::string(80,'#') << std::endl;
  std::cout << "### Integration on the unit interval [0,1] ###" << std::endl;
  std::cout << std::endl;
  std::cout << "Gauss-Legendre quadrature of order "<< 2*N << std::endl;
  std::cout << std::endl;
  std::cout << "Test monomial P(x) = x**p" << std::endl;
  std::cout << std::endl;
  for (int p = 0; p <= 2*N; p++) {
    auto Pp = [p] (Eigen::VectorXd x) {return pow(x(0),p);};
    double Q = integrate(QR1D,Pp);
    double I = 1./(p + 1.);
    double Err = std::abs(I - Q);
    std::cout << "  Absolute error (p = " << p << ") : "
              << Err << std::endl;
  }

  // test 2D quadrilateral quadrature on monomials in unit quadrilateral
  std::cout << std::string(80,'#') << std::endl;
  std::cout << "### Integration on the unit quadrilateral [0,1]**2 ###"
            << std::endl;
  std::cout << std::endl;
  std::cout << "2D quadrature of order "<< 2*N << std::endl;
  std::cout << std::endl;
  std::cout << "Test monomial P(x1,x2) = x1**p * x2**q" << std::endl;
  std::cout << std::endl;
  for (int i = 0; i <= 2*N; i++) {
    for (int j = 0; j <= i; j++) {
      int p = i - j;
      int q = j;
      auto Ppq = [p,q] (Eigen::VectorXd x) {return pow(x(0),p)*pow(x(1),q);};
      double Q = integrate(Quad_QR2D,Ppq);
      double I = 1./((p + 1.)*(q + 1.));
      double Err = std::abs(I - Q);
      std::cout << "  Absolute error (p = " << p << ", q = " << q
                << ", p + q = " << p + q << ") : "
                << Err << std::endl;
    }
  }

  // test 2D triangular quadrature on monomials in unit triangle
  std::cout << std::string(80,'#') << std::endl;
  std::cout << "### Integration on the unit triangle K ###"
            << std::endl;
  std::cout << std::endl;
  std::cout << "2D quadrature of order "<< 2*N - 1 << std::endl;
  std::cout << std::endl;
  std::cout << "Test monomial P(x1,x2) = x1**p * x2**q" << std::endl;
  std::cout << std::endl;
  for (int i = 0; i <= 2*N; i++) {
    for (int j = 0; j <= i; j++) {
      int p = i - j;
      int q = j;
      auto Ppq = [p,q] (Eigen::VectorXd x) {return pow(x(0),p)*pow(x(1),q);};
      double Q = integrate(Tria_QR2D,Ppq);
      double I = double(factorial(p)*factorial(q))
               / double(factorial(p + q + 2));
      double Err = std::abs(I - Q);
      std::cout << "  Absolute error (p = " << p << ", q = " << q
                << ", p + q = " << p + q << ") : "
                << Err << std::endl;
    }
  }

  return 0;

}
