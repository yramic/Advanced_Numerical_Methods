/**
 * @file tensorproductchebintp_test.cc
 * @brief ADVNCSE homework TensorProductChebIntp test code
 * @author Bob Schreiner
 * @date August 2023
 * @copyright Developed at SAM, ETH Zurich
 */

#include "../tensorproductchebintp.h"

#include <gtest/gtest.h>

namespace TensorProductChebIntp::test {
TEST(TensorProductChebIntp, chebInterpEval1D) {
  auto kernel = [](double x) {
    return 1. / (1. + (x * x));
  };  // the kernel we want to interpolate

  // We create a random number generator to sample points between [-1 , 1]
  std::random_device rand_dev;
  std::mt19937 generator(rand_dev());
  std::uniform_real_distribution<double> distr(-1., 1.);

  // We sample points between [-1 , 1]
  const int N = 100;
  std::vector<double> x(N);
  for (int n = 0; n < N; ++n) x[n] = distr(generator);

  // Compute the error for the degree 2^q
  double error;
  const int q = 16;
  const std::vector<double> res =
      TensorProductChebIntp::chebInterpEval1D(q, kernel, x);
  error = TensorProductChebIntp::errorestimate(q, kernel, x, res);
  ASSERT_NEAR(error, 0.0, 1e-8);
}

TEST(TensorProductChebIntp, chebInterpEval2D) {
  auto kernel = [](double x, double y) {
    return 1. / (1. + (x - y) * (x - y));
  };  // the kernel we want to interpolate

  // We create a random number generator to sample points between [-1 , 1]
  std::random_device rand_dev;
  std::mt19937 generator(rand_dev());
  std::uniform_real_distribution<double> distr(-1., 1.);

  // We sample points between [-1 , 1]
  const int N = 100;
  std::vector<Eigen::Vector2d> x(N);
  for (int n = 0; n < N; ++n)
    x[n] = (Eigen::Vector2d() << distr(generator), distr(generator)).finished();

  // Compute the error for the degree 2^q
  double error;
  const int q = 16;
  const std::vector<double> res =
      TensorProductChebIntp::chebInterpEval2D(q, kernel, x);
  error = TensorProductChebIntp::errorestimate(q, kernel, x, res);
  ASSERT_NEAR(error, 0.0, 1e-8);
}

TEST(TensorProductChebIntp, chebInterpEval2D_exact) {
  auto kernel = [](double x, double y) {
    return 1. / (1. + (x - y) * (x - y));
  };  // the kernel we want to interpolate

  // We create a random number generator to sample points between [-1 , 1]
  std::random_device rand_dev;
  std::mt19937 generator(rand_dev());
  std::uniform_real_distribution<double> distr(-1., 1.);

  // We sample points between [-1 , 1]
  const int q = 16;
  const int N = q * q;
  std::vector<Eigen::Vector2d> x(N);
  for (int n = 0; n < q; ++n)
    for (int k = 0; k < q; ++k)
      x[n] = (Eigen::Vector2d()
                  << std::cos((2.0 * (n + 1) - 1.0) / (2 * q) * M_PI),
              std::cos((2.0 * (k + 1) - 1.0) / (2 * q) * M_PI))
                 .finished();

  // Compute the error for the degree 2^q
  double error;

  const std::vector<double> res =
      TensorProductChebIntp::chebInterpEval2D(q, kernel, x);
  error = TensorProductChebIntp::errorestimate(q, kernel, x, res);
  ASSERT_NEAR(error, 0.0, 1e-10);
}

TEST(TensorProductChebIntp, genChebInterpEval2D) {
  auto kernel = [](double x, double y) {
    return 1. / (1. + (x - y) * (x - y));
  };  // the kernel we want to interpolate

  // We create a random number generator to sample points between [-1 , 1]
  std::random_device rand_dev;
  std::mt19937 generator(rand_dev());
  Eigen::Vector2d a;
  Eigen::Vector2d b;
  const int N = 100;

  a << -1.6, 0.5;
  b << -2.5, -0.6;
  std::uniform_real_distribution<double> distr2(a[0], b[0]);
  std::uniform_real_distribution<double> distr3(a[1], b[1]);
  std::vector<Eigen::Vector2d> x2(N);
  for (int n = 0; n < N; ++n)
    x2[n] =
        (Eigen::Vector2d() << distr2(generator), distr3(generator)).finished();

  double error;
  int q = 16;
  const std::vector<double> res2 =
      TensorProductChebIntp::genChebInterpEval2D(q, kernel, a, b, x2);
  error = TensorProductChebIntp::errorestimate(q, kernel, x2, res2);
  ASSERT_NEAR(error, 0.0, 1e-8);
}
}  // namespace TensorProductChebIntp::test
