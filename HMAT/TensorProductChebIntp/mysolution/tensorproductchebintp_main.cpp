/**
 * @ file tensorproductchebintp_main.cpp
 * @ brief ADVNCSE homework TensorProductChebIntp MAIN FILE
 * @ author Bob Schreiner
 * @ date August 2023
 * @ copyright Developed at SAM, ETH Zurich
 */

#include "tensorproductchebintp.h"
#include <random>

int main(int /*argc*/, char** /*argv*/) {
  auto kernel1D = [](double x){return 1./(1. + (x*x));}; // the 1d kernel we want to interpolate
  auto kernel2D = [](double x , double y){return 1./(1. + (x-y)*(x-y));}; // the 2d kernel we want to interpolate

  // We create a random number generator to sample points between [-1 , 1]
  std::random_device                  rand_dev;
  std::mt19937                        generator(rand_dev());
  std::uniform_real_distribution<double>  distr(-1., 1.);

    // We sample points between [-1 , 1]^2
    const int N = 100;
    std::vector<double> y(N);
    for (int n=0; n<N; ++n)
        y[n] = distr(generator);

    // We tabulate the error
    double error;
    std::cout << "Testing k(x,y) = 1./(1. + (x*x))" << std::endl;
    std::cout << "q" << std::setw(20) << "error" << std::endl;
    for (int q = 1; q<10;q++){
        const std::vector<double> res = TensorProductChebIntp::chebInterpEval1D(std::pow(2,q), kernel1D, y);
        error = TensorProductChebIntp::errorestimate(std::pow(2,q), kernel1D, y, res);
        std::cout << std::pow(2,q) << std::setw(20) << error << std::endl;
    }

  // We sample points between [-1 , 1]^2
  std::vector<Eigen::Vector2d> x(N);
  for (int n=0; n<N; ++n)
    x[n] = (Eigen::Vector2d() << distr(generator), distr(generator)).finished();

  // We tabulate the error
  std::cout << "\n\nTesting k(x,y) = 1./(1. + (x-y)*(x-y))" << std::endl;
  std::cout << "q" << std::setw(20) << "error" << std::endl;
  for (int q = 1; q<10;q++){
    const std::vector<double> res = TensorProductChebIntp::chebInterpEval2D(std::pow(2,q), kernel2D, x);
    error = TensorProductChebIntp::errorestimate(std::pow(2,q), kernel2D, x, res);
    std::cout << std::pow(2,q) << std::setw(20) << error << std::endl;
  }

  // We sample points in the real numbers
  Eigen::Vector2d a;
  Eigen::Vector2d b;

  a << -1.6 , 0.5;
  b << -2.5 , -0.6;
  std::uniform_real_distribution<double>  distr2(a[0] , b[0]);
  std::uniform_real_distribution<double>  distr3(a[1] , b[1]);
  std::vector<Eigen::Vector2d> x2(N);
  for (int n=0; n<N; ++n)
    x2[n] = (Eigen::Vector2d() << distr2(generator), distr3(generator)).finished();

  // We tabulate the error
  std::cout << "\n\nTesting k(x,y) = 1./(1. + (x-y)*(x-y))" << std::endl;
  std::cout << "q" << std::setw(20) << "error" << std::endl;
  for (int q = 1; q<10;q++){
    const std::vector<double> res2 = TensorProductChebIntp::genChebInterpEval2D(std::pow(2,q), kernel2D,a,b, x2);
    error = TensorProductChebIntp::errorestimate(std::pow(2,q), kernel2D, x2 , res2);
    std::cout << std::pow(2,q) << std::setw(20) << error << std::endl;
  }
  return 0;
}
