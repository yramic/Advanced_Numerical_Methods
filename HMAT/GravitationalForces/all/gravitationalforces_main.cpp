/**
 * @ file gravitationalforces_main.cpp
 * @ brief NPDE homework XXX MAIN FILE
 * @ author
 * @ date
 * @ copyright Developed at SAM, ETH Zurich
 */

#include <string_view>
#include <vector>

#include "gravitationalforces.h"

using namespace std::literals;

int main(int argc, char **argv) {
  std::cout << "GravitationalForces ADVNCSE HW projects" << std::endl;
  // Number of stars can be given through the command line
  unsigned int n_def = 1000;  // Default
  if (argc > 2) {
    if (argv[1] == "-n"sv || argv[1] == "-s"sv || argv[1] == "-stars"sv) {
      n_def = strtol(argv[2], nullptr, 0);
    }
  }

  /*
  {
    // Simple test case for the computation of exact forces
    std::vector<Eigen::Vector2d> pos(4, Eigen::Vector2d::Zero());
    const std::vector<double> mass(4, 1.0);
    pos[0] << 0.0, 0.0;
    pos[1] << 0.0, 1.0;
    pos[2] << 1.0, 0.0;
    pos[3] << 1.0, 1.0;
    std::vector<Eigen::Vector2d> forces =
        GravitationalForces::computeForces_direct(pos, mass);

    std::cout << "Forces:" << std::endl;
    for (Eigen::Vector2d &f : forces) {
      std::cout << "[" << f.transpose() << "] ";
    }
    std::cout << std::endl;
  }
  */
  /*
  {
    // "Initialize stars"
    std::vector<Eigen::Vector2d> pos =
        GravitationalForces::initStarPositions(n_def);
    const std::vector<double> mass(n_def, 1.0);
    // Initialize quadtree of star clusters
    const GravitationalForces::StarQuadTreeClustering clustering(pos, mass);
    // clustering.outputQuadTree(std::cout);

    // Admissibility parameters to be investigated
    const std::vector<double> etas = {1.0,  1.25, 1.5,  1.75, 2.0,
                                      2.25, 2.5,  2.75, 3.0};
    // Compute force errors
    std::vector<double> errs =
        GravitationalForces::forceError(clustering, etas);
    std::cout << std::setw(16) << "eta" << std::setw(16) << "error"
              << std::endl;
    for (int l = 0; l < etas.size(); ++l) {
      std::cout << std::setw(16) << etas[l] << std::setw(16) << errs[l]
                << std::endl;
    }
  }
  */
  // Measure runtimes
  unsigned int n = 20;
  std::cout << std::setw(10) << "n" << std::setw(10) << "runtime(exact)"
            << std::setw(10) << "runtime(clustering)" << std::endl;
  for (int l = 0; l < 10; n *= 2, ++l) {
    auto [ms_exact, ms_cluster] = GravitationalForces::measureRuntimes(n);
    std::cout << std::setw(10) << n << std::setw(10) << ms_exact
              << std::setw(10) << ms_cluster << std::endl;
  }
  return 0;
}
