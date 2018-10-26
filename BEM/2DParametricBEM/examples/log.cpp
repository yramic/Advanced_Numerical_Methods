#include "integral_gauss.hpp"
#include <Eigen/Dense>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <fstream>
#include <string>

int main() {
  unsigned int num_pts = 100;
  std::string filename = "log_out.txt";
  std::ofstream output(filename);
  output << std::setw(15) << "#y" << std::setw(15) << "exact" << std::setw(15)
         << "calculated" << std::setw(15) << "order" << std::endl;
  Eigen::VectorXd y1 = Eigen::VectorXd::LinSpaced(num_pts/2, -2.5, -1.01);
  Eigen::VectorXd y2 = Eigen::VectorXd::LinSpaced(num_pts/2, -1.01, -1.001);
  Eigen::VectorXd y(num_pts);
  y << y1,y2;
  Eigen::VectorXd analytic_values(num_pts);
  for (unsigned int i = 0; i < num_pts; ++i)
    analytic_values(i) =
        (1 - y(i)) * (log(1 - y(i)) - 1) + (1 + y(i)) * (log(-1 - y(i)) - 1);

  Eigen::VectorXd computed_values(num_pts);
  Eigen::VectorXd convergence_order(num_pts);
  for (unsigned int i = 0; i < num_pts; ++i) {
    std::function<double(double)> integrand = [&](double x) {
      return log(x - y(i));
    };
    // Find the converged integral and convergence order
    int order = -3;
    double integral_old, integral_new = 0.;
    do {
      order += 5;
      integral_old = integral_new;
      integral_new =
          parametricbem2d::ComputeIntegral(integrand, -1., 1., order);
    } while (fabs(integral_new - integral_old) > 1e-9 &&
             order < 1000); // Convergence criteria
    if (order > 1000)
      std::cout << "Not converged for y = " << y(i) << std::endl;

    computed_values(i) = integral_new;
    convergence_order(i) = order;
  }
  for (unsigned int i = 0; i < num_pts; ++i)
    output << std::setw(15) << y(i) << std::setw(15) << analytic_values(i)
           << std::setw(15) << computed_values(i) << std::setw(15)
           << convergence_order(i) << std::endl;
  return 0;
}
