/**
 * @ file gravitationalforces_main.cpp
 * @ brief NPDE homework XXX MAIN FILE
 * @ author
 * @ date
 * @ copyright Developed at SAM, ETH Zurich
 */

#include "gravitationalforces.h"

int main(int /*argc*/, char** /*argv*/) {
  // Initialize 100 stars
  const int n = 100;
  std::vector<Eigen::Vector2d> pos =
      GravitationalForces::initStarPositions(n);
  // Output star positions
  std::cout << "Star positions:" << std::endl;
  for (Eigen::Vector2d &p : pos) {
    std::cout << "[" << p.transpose() << "] ";
  }
  std::cout << std::endl;
  
  runtimemeasuredemo();
  
  return 0;
}
