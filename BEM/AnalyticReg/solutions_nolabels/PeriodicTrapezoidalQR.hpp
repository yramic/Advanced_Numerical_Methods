////
//// Copyright (C) 2017 SAM (D-MATH) @ ETH Zurich
//// Author(s): curzuato < >
//// Contributors:  dcasati
//// This file is part of the AdvNumCSE repository.
////
#ifndef PTRAPEZOIDALQR_HPP
#define PTRAPEZOIDALQR_HPP

#include <Eigen/Dense>
#include <cmath>
#include <iostream>

/* @brief Compute Periodic Trapezoidal Rule over I=[0,2 PI]
 * \param[in] n Number of quadrature points
 */
std::pair<Eigen::VectorXd, double> PeriodicTrapRule(int N) {
  double wq = 2. * M_PI / N;

  Eigen::VectorXd aux;
  aux.setLinSpaced(N + 1, 0, 1);
  Eigen::VectorXd xq = aux.segment(0, N) * 2. * M_PI;

  // return
  return std::make_pair(xq, wq);
}

#endif
