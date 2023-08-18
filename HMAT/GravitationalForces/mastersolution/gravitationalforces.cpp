/**
 * @file GravitationalForces.cpp
 * @brief NPDE homework GravitationalForces code
 * @author
 * @date
 * @copyright Developed at SAM, ETH Zurich
 */

#include "gravitationalforces.h"

#include <chrono>

namespace GravitationalForces {

/* SAM_LISTING_BEGIN_1 */
std::vector<Eigen::Vector2d> initStarPositions(unsigned int n, double mindist) {
  assertm(mindist > 0, "Minimal distance must be positive");
  const double max_md = 1.0 / (2.0 * std::sqrt(n));  // Maximal minimal distance
  if (mindist > max_md) {
    mindist = max_md;
  }
  // Vector for returning positions
  std::vector<Eigen::Vector2d> pos;
  // Fill position vector with random positions taking into account minimal
  // distance requirement
  do {
    Eigen::Matrix<double, 2, Eigen::Dynamic> randpos =
        0.5 * (Eigen::Matrix<double, 2, Eigen::Dynamic>::Random(2, n) +
               Eigen::Matrix<double, 2, Eigen::Dynamic>::Constant(2, n, 1.0));
    // Add columns of random matrix one by one checking minimal distance
    for (int i = 0; (i < n) && (pos.size() < n); ++i) {
      int m = 0;
      bool good = true;
      while ((m < pos.size()) and
             (good = (pos[m] - randpos.col(i)).norm() > mindist)) {
        ++m;
      }
      if (good) {
        // Accept next position vector
        pos.emplace_back(randpos.col(i));
      }
      // It can happen that positions are rejected. In this case we have to try
      // more.
    }
  } while (pos.size() < n);
  return pos;
}
/* SAM_LISTING_END_1 */

/* SAM_LISTING_BEGIN_2 */
std::vector<Eigen::Vector2d> computeForces_direct(
    const std::vector<Eigen::Vector2d> &masspositions,
    const std::vector<double> &masses) {
  const unsigned int n = masspositions.size();
  assertm(n == masses.size(),
          "Mismathc of sizes of masspositions and masses vectors");
  std::vector<Eigen::Vector2d> forces(n);

  return forces;
}
/* SAM_LISTING_END_2 */

StarQuadTree::StarQuadTree(const std::vector<Eigen::Vector2d> &starpos,
                           const std::vector<double> &starmasses)
    : n(starpos.size()), starpos_(starpos), starmasses_(starmasses) {
  std::vector<unsigned int> star_idx(n);
  for (int i = 0; i < n; ++i) {
    star_idx[i] = i;
  }
  root_ = std::move(
      std::make_unique<StarQuadTree::StarQuadTreeNode>(star_idx, *this));
}

StarQuadTree::StarQuadTreeNode::StarQuadTreeNode(
    std::vector<unsigned int> star_idx, const StarQuadTree &tree)
    : star_idx_(star_idx) {}

StarQuadTreeClustering::StarQuadTreeClustering(
    const std::vector<Eigen::Vector2d> &starpos,
    const std::vector<double> &starmasses)
    : StarQuadTree(starpos, starmasses) {}

/* SAM_LISTING_BEGIN_5 */
bool StarQuadTreeClustering::isAdmissible(const StarQuadTreeNode &node,
                                          Eigen::Vector2d p, double eta) const {
  return true;
}
/* SAM_LISTING_END_5 */

/* SAM_LISTING_BEGIN_6 */
Eigen::Vector2d StarQuadTreeClustering::forceOnStar(unsigned int j,
                                                    double eta) const {}
/* SAM_LISTING_END_6 */

/* SAM_LISTING_BEGIN_7 */
std::vector<double> forceError(const StarQuadTreeClustering &qt,
                               const std::vector<double> &etas);
/* SAM_LISTING_END_7 */

}  // namespace GravitationalForces

/* SAM_LISTING_BEGIN_8 */
void runtimemeasuredemo(void) {
  auto t1 = std::chrono::high_resolution_clock::now();
  double s = 0.0;
  for (long int i = 0; i < 10000000; ++i) {
    s += 1.0 / std::sqrt((double)i);
  }
  auto t2 = std::chrono::high_resolution_clock::now();
  /* Getting number of milliseconds as a double. */
  std::chrono::duration<double, std::milli> ms_double = t2 - t1;
  std::cout << "Runtime = " << ms_double.count() << "ms\n";
}
/* SAM_LISTING_END_8 */
