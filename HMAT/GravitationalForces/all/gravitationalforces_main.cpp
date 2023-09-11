/**
 * @ file gravitationalforces_main.cpp
 * @ brief NPDE homework XXX MAIN FILE
 * @ author
 * @ date
 * @ copyright Developed at SAM, ETH Zurich
 */

#include <string_view>

#include "gravitationalforces.h"

using namespace GravitationalForces;
using namespace std::literals;

int main(int argc, char **argv) {
  // default number of stars (can be adjusted using args "-n", "-s" or,
  // "-stars")
  unsigned int n = 320;

  // default mode (can be adjusted using args "-t", "-timing" or "-e", "-error")
  unsigned int timing_mode = 0;

  if (argc > 2) {
    for (int i = 1; i < argc; i += 2) {
      if (argv[i] == "-n"sv || argv[i] == "-s"sv || argv[i] == "-stars"sv) {
        n = strtol(argv[i + 1], nullptr, 0);
      }
      if (argv[i] == "-t"sv || argv[i] == "-e"sv || argv[i] == "-timing"sv ||
          argv[i] == "-error"sv) {
        timing_mode = strtol(argv[i + 1], nullptr, 0);
      }
    }
  }

  // Initialize stars
  std::vector<Eigen::Vector2d> pos = GravitationalForces::initStarPositions(n);
  const std::vector<double> mass(n, 1.0);

  /* Output star positions
  std::cout << "Star positions:" << std::endl;
  for (Eigen::Vector2d &p : pos) {
    std::cout << "[" << p.transpose() << "] ";
  }
  std::cout << std::endl;
  */

  /* Test forces
  std::vector<Eigen::Vector2d> pos_2(n, Eigen::Vector2d::Zero());
  pos_2[0] << 0.0, 0.0;
  pos_2[1] << 0.0, 1.0;
  pos_2[2] << 1.0, 0.0;
  pos_2[3] << 1.0, 1.0;
  std::vector<Eigen::Vector2d> forces =
      GravitationalForces::computeForces_direct(pos_2, mass);

  std::cout << "Forces:" << std::endl;
  for (Eigen::Vector2d &f : forces) {
    std::cout << "[" << f.transpose() << "] ";
  }
  std::cout << std::endl;
  */
  // Test QuadTree
  StarQuadTree *const tree = new StarQuadTree(pos, mass);
  std::cout << std::endl << "All Indices: ";
  for (int i = 0; i < tree->root_->star_idx_.size(); i++) {
    if (i == 10) std::cout << "... , ";
    if (i < 10 || i > tree->root_->star_idx_.size() - 11)
      std::cout << tree->root_->star_idx_[i] << ", ";
  }
  std::cout << std::endl;

  for (unsigned int i : {0, 1, 2, 3}) {
    std::cout << std::endl << "SON " << i << std::endl;
    if (!tree->root_->sons_[i]) continue;
    std::cout << "Indices: ";
    for (int j = 0; j < tree->root_->sons_[i]->star_idx_.size(); j++) {
      if (j == 10) std::cout << "... , ";
      if (j < 10 || j > tree->root_->sons_[i]->star_idx_.size() - 11)
        std::cout << tree->root_->sons_[i]->star_idx_[j] << ", ";
    }
    std::cout << std::endl;
    for (unsigned int k : {0, 1, 2, 3}) {
      std::cout << "\t SON " << i << k << ":  ";
      if (tree->root_->sons_[i]->sons_[k] == 0) {
        std::cout << std::endl << std::endl;
        continue;
      }
      for (int ik = 0; ik < tree->root_->sons_[i]->sons_[k]->star_idx_.size();
           ik++) {
        if (ik == 5)
          std::cout << "... " << std::endl
                    << "\t\t  ... " << std::endl
                    << "\t\t  ";
        if (ik < 5 ||
            ik > tree->root_->sons_[i]->sons_[k]->star_idx_.size() - 6)
          std::cout
              << "Idx: " << tree->root_->sons_[i]->sons_[k]->star_idx_[ik]
              << "\t POS: "
              << tree->starpos_[tree->root_->sons_[i]->sons_[k]->star_idx_[ik]]
                               [0]
              << ", "
              << tree->starpos_[tree->root_->sons_[i]->sons_[k]->star_idx_[ik]]
                               [1]
              << std::endl
              << "\t\t  ";
      }
      std::cout << "MASS: " << tree->root_->sons_[i]->sons_[k]->mass
                << "\t CENTER: " << tree->root_->sons_[i]->sons_[k]->center[0]
                << ", " << tree->root_->sons_[i]->sons_[k]->center[1]
                << std::endl
                << std::endl;
    }
  }

  StarQuadTreeClustering *const clustering =
      new StarQuadTreeClustering(pos, mass);

  std::vector<double> etas;  // = {0.05,0.06,0.07,0.08,0.09,0.1};

  // Work with multiple etas
  for (double i = 0.8; i <= 2.05; i += 0.1) {
    etas.push_back(i);
  }

  // Work with a single eta or comment out the following single line
  // etas.clear(); etas.push_back(1);

  etas = measureRuntimes(*clustering, etas, timing_mode);

  std::cout << std::endl << "Errors(eta): ";
  std::stringstream ss;

  for (auto it = etas.begin(); it != etas.end(); it++) {
    if (it != etas.begin()) {
      ss << ", ";
    }
    ss << *it;
  }
  std::cout << ss.str() << std::endl;

  // runtimemeasuredemo();

  return 0;
}
