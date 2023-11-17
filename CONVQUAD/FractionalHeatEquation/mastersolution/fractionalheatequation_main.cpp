/**
 * @ file fractionalheatequation_main.cpp
 * @ brief NPDE homework FractionalHeatEquation MAIN FILE
 * @ author Jörg Nick, Bob Schreiner
 * @ date October 2023
 * @ copyright Developed at SAM, ETH Zurich
 */

#include <chrono>
#include <fstream>

#include "fractionalheatequation.h"

/* SAM_LISTING_BEGIN_0 */
int main(int /*argc*/, char** /*argv*/) {
  FractionalHeatEquation::SqrtsMplusA Amat(2, std::complex<double>(1, 1));
  std::cout << "Running code for ADVNCSE HW FractionalHeatEquation"
            << std::endl;
  const double T = 1.0;
  const int L_start = 4;
  const int L_max = 13;
  const int n = 4;
  std::function<double(double, Eigen::Vector2d)> f =
      [](double t, Eigen::Vector2d x) { return t * t * t; };

  const unsigned num_repetitions = 5;
  std::ofstream out(CURRENT_SOURCE_DIR "/runtimes.csv");

  std::cout << std::left << std::setw(10) << "M" << std::left << std::setw(15)
            << "time_MOT" << std::left << std::setw(15) << "time_TOEP"
            << std::left << std::setw(15) << "time_ASAO" << std::left
            << std::setw(15) << std::endl;
  out << "M"
      << ","
      << "time_MOT"
      << ","
      << "time_TOEP"
      << ","
      << "time_ASAO" << std::endl;
  for (unsigned int L = L_start; L < L_max; ++L) {
    double time_mot = std::numeric_limits<double>::max();
    double time_toep = std::numeric_limits<double>::max();
    double time_asao = std::numeric_limits<double>::max();
    for (int k = 0; k < num_repetitions; k++) {
      // Use C++ chrono library to measure runtimes
      auto t1 = std::chrono::high_resolution_clock::now();
      Eigen::VectorXd mu_MOT =
          FractionalHeatEquation::evlMOT(f, n, T, std::pow(2, L) - 1);
      auto t2 = std::chrono::high_resolution_clock::now();
      // Getting number of seconds as a double
      std::chrono::duration<double> ms_mot = (t2 - t1);

      time_mot = std::min(time_mot, ms_mot.count());

      t1 = std::chrono::high_resolution_clock::now();
      Eigen::VectorXd mu_Toep =
          FractionalHeatEquation::evlTriangToeplitz(f, n, T, L);
      t2 = std::chrono::high_resolution_clock::now();
      std::chrono::duration<double> ms_toep = (t2 - t1);

      // Taking the minimal measured time as the result
      time_toep = std::min(time_toep, ms_toep.count());

      t1 = std::chrono::high_resolution_clock::now();
      Eigen::VectorXd mu_ASAO = FractionalHeatEquation::evlASAOCQ(f, n, T, L);

      t2 = std::chrono::high_resolution_clock::now();
      // Getting number of seconds as a double
      std::chrono::duration<double> ms_asao = (t2 - t1);

      time_asao = std::min(time_asao, ms_asao.count());
    }
    std::cout << std::left << std::setw(10) << std::pow(2, L) - 1 << std::left
              << std::setw(15) << time_mot << std::left << std::setw(15)
              << time_toep << std::left << std::setw(15) << time_asao
              << std::left << std::setw(15) << std::endl;
    out << std::pow(2, L) - 1 << "," << time_mot << "," << time_toep << ","
        << time_asao << std::endl;
  }
  out.close();
  // Call python script
  std::system("python3 " CURRENT_SOURCE_DIR "/plot.py " CURRENT_SOURCE_DIR
              "/runtimes.csv " CURRENT_SOURCE_DIR "/runtimes.png");
  return 0;
}
/* SAM_LISTING_END_0 */