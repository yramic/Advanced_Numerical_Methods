/**
 * @ file stationarylineariterations_main.cpp
 * @ brief ACVNCSE homework StationaryLinearIterations MAIN FILE
 * @  Bob Schreiner
 * @ date
 * @ copyright Developed at SAM, ETH Zurich
 */

#include <fstream>
#include <iomanip>

#include "stationarylineariterations.h"

int main(int /*argc*/, char** /*argv*/) {
  std::cout << "Running code for ADVNCSE HW StationaryLinearIterations"
            << std::endl;
  const unsigned c_max = 11;
  const unsigned l_max = 11;
  Eigen::MatrixXd results(c_max, l_max);
  std::cout << "Asymptotic rate of convergence of: Gauss Seidel" << std::endl;
  std::cout << std::left << std::setw(10) << "\nc " << std::left
            << std::setw(10) << "n" << std::left << std::setw(10) << "lambda(X)"
            << std::endl;

  std::ofstream out(CURRENT_SOURCE_DIR "/convergence.csv");
  out << "c,n,lambda(X)" << std::endl;
  for (unsigned c = 0; c < c_max; ++c) {
    for (unsigned l = 1; l < l_max; ++l) {
      results(c, l) =
          StationaryLinearIterations::gaussSeidelRate(std::pow(2, l), c);
      std::cout << std::left << std::setw(10) << c << std::left << std::setw(10)
                << std::pow(2, l) << std::left << std::setw(10) << results(c, l)
                << std::endl;
      out << c << "," << std::pow(2, l) << "," << results(c, l) << std::endl;
    }
  }
  out.close();
  // Call python script
  std::system("python3 " CURRENT_SOURCE_DIR "/plot.py " CURRENT_SOURCE_DIR
              "/convergence.csv " CURRENT_SOURCE_DIR "/convergence.png");
  return 0;
}
