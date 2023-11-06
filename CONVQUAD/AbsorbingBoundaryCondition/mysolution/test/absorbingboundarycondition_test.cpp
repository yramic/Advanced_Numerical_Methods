/**
 * @file absorbingboundarycondition_test.cc
 * @brief ADVNCSE homework AbsorbingBoundaryCondition test code
 * @author Peiyuan Xie
 * @date November 2023
 * @copyright Developed at SAM, ETH Zurich
 */

#include "../absorbingboundarycondition.h"

#include <gtest/gtest.h>

namespace AbsorbingBoundaryCondition::test {

// Test cqweights_by_dft()
TEST(AbsorbingBoundaryCondition, test_cqweights) {
  // For $\delta(z)=1-z, F(s)=s$,
  // exact solution: $CQ_\tau(F)=(1/\tau, -1/\tau, 0, ..., 0)\in R^{M+1}$
  auto F = [](std::complex<double> s) { return s; };
  auto delta = [](std::complex<double> z) {
    return std::complex<double>(1) - z;
  };
  double tau = 1;
  size_t N = 10;
  Eigen::VectorXd w = cqweights_by_dft(F, delta, tau, N);
  // Gound truth
  Eigen::VectorXd gt(w.size());
  gt.setZero();
  gt(0) = 1;
  gt(1) = -1;

  ASSERT_NEAR((gt - w).norm(), 0.0, 1e-8);
}

// Test solve_IBVP()
TEST(AbsorbingBoundaryCondition, test_solve) {
  auto g = [](double t) { return sin(M_PI * t); };
  // Implemented solution
  int M_ref = 512;
  int N_ref = 512;
  double T0 = 0.5;
  VectorXd u0 = solve_IBVP(g, M_ref, N_ref, T0);
  double T1 = 1;
  VectorXd u1 = solve_IBVP(g, M_ref, N_ref, T1);
  // Recover Neumann boundary condition at x=0
  ASSERT_NEAR(N_ref * (u0(1) - u0(0)), 1.0, 1e-5);
  ASSERT_NEAR(N_ref * (u1(1) - u1(0)), 0.0, 1e-5);
}

}  // namespace AbsorbingBoundaryCondition::test
