#include "../abelintegralequation.h"

#include <gtest/gtest.h>
namespace AbelIntegralEquation::test {
// Test case
auto u = [](double t) { return 2. / M_PI * sqrt(t); };
auto y = [](double t) { return t; };

TEST(AbelIntegralEquation, test_poly_spec) {
  int p = 10;
  double tau = 0.1;
  size_t N = round(1. / tau);
  VectorXd u_ex(N + 1);
  // Reference solution
  u_ex << 0.0284021, 0.201482, 0.284706, 0.348974, 0.402129, 0.450595, 0.493023,
      0.53237, 0.569899, 0.603342, 0.631822;
  // Implemented solution
  VectorXd u_app = AbelIntegralEquation::poly_spec_abel(y, p, tau);
  // Difference between reference and implementation
  double error = (u_app - u_ex).norm();

  ASSERT_NEAR(error, 0.0, 1e-5);
}

TEST(AbelIntegralEquation, test_cq_ieul) {
  size_t N = 10;
  VectorXd u_ex(N + 1);
  // Reference solution
  u_ex << 0, 0.178412, 0.267619, 0.334523, 0.390277, 0.439062, 0.482968,
      0.523215, 0.560588, 0.595625, 0.628715;
  // Implemented solution
  VectorXd u_app = AbelIntegralEquation::cq_ieul_abel(y, N);
  // Difference between reference and implementation
  double error = (u_app - u_ex).norm();

  ASSERT_NEAR(error, 0.0, 1e-5);
}

TEST(AbelIntegralEquation, test_cq_bdf2) {
  size_t N = 10;
  VectorXd u_ex(N + 1);
  // Reference solution
  u_ex << 0, 0.21851, 0.291346, 0.352043, 0.404648, 0.451519, 0.49412, 0.533404,
      0.570029, 0.604462, 0.637052;
  // Implemented solution
  VectorXd u_app = AbelIntegralEquation::cq_bdf2_abel(y, N);
  // Difference between reference and implementation
  double error = (u_app - u_ex).norm();

  ASSERT_NEAR(error, 0.0, 1e-5);
}

}  // namespace AbelIntegralEquation::test