/**
 * @ file gravitationalforces_main.cpp
 * @ brief NPDE homework GravitationalForces MAIN FILE
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

  {
    // "Initialize stars"
    std::vector<Eigen::Vector2d> pos =
        GravitationalForces::initStarPositions(n_def);
    const std::vector<double> mass(n_def, 1.0);
    // Initialize quadtree of star clusters
    const GravitationalForces::StarQuadTreeClustering clustering(pos, mass);
    // clustering.outputQuadTree(std::cout);

    //--------------------------------------------------------------------------
    // This was implemented to Test my Solution on the forces for 9 Iterations:
    // std::vector<Eigen::Vector2d> forces =
    //     GravitationalForces::computeForces_direct(pos, mass);
    
    // for(unsigned int i {0}; i < 10; ++i) {
    //   std::cout << "Force from " << i << ": " << forces[i] << std::endl;
    // }
    //--------------------------------------------------------------------------

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
