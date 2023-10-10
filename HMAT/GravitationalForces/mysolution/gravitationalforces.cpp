/**
 * @file GravitationalForces.cpp
 * @brief NPDE homework GravitationalForces code
 * @author
 * @date
 * @copyright Developed at SAM, ETH Zurich
 */

#include "gravitationalforces.h"

#include <Eigen/src/Core/Matrix.h>

#include <chrono>
#include <limits>

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
// **********************************************************************

// int counter {1}; // I want to count how many Iterations the code takes to reach a solution

do {
  // First I want to generate a 2xn Matrix with random values in between [0, 1]:
  Eigen::Matrix<double, 2, Eigen::Dynamic> rand_pos {
    0.5* (Eigen::Matrix<double, 2, Eigen::Dynamic>::Random(2,n) +
          Eigen::Matrix<double, 2, Eigen::Dynamic>::Constant(2, n, 1.0))
    };
    for (int i {0}; (i < n) && (pos.size() < n); ++i) {
      int j {0};
      bool test {true};
      while ((j < pos.size()) and (test = (pos[j] - rand_pos.col(i)).norm() > mindist)) {
        ++j;
      }
      if (test) {
        pos.emplace_back(rand_pos.col(i));
      }
    }
  } while(pos.size() < n);

  // The following code snipped was implemented just for testing!
  // for (int i {0}; i < 10; ++i) {
  //   std::cout<< "Element " << i << " :\n" << pos[i] << std::endl; 
  // }

  // This gives back the number of Iterations necessary to reach a solution:
  // std::cout<< "Iteration: " << counter << std::endl;
  // counter++;

// **********************************************************************
  return pos;
}
/* SAM_LISTING_END_1 */

/* SAM_LISTING_BEGIN_2 */
std::vector<Eigen::Vector2d> computeForces_direct(
    const std::vector<Eigen::Vector2d> &masspositions,
    const std::vector<double> &masses) {
  const unsigned int n = masspositions.size();
  assertm(n == masses.size(),
          "Mismatch of sizes of masspositions and masses vectors");
  // Forces will be stored in this array
  std::vector<Eigen::Vector2d> forces(n);
// **********************************************************************
  Eigen::Vector2d acc;
  for (unsigned int j = 0; j < n; j++) {
    // Compute force on $\cob{j}$-th star
    acc.setZero();
    for (unsigned int i = 0; i < n; i++) {
      if (i != j) {
        const Eigen::Vector2d diff_masspositions =
            masspositions[i] - masspositions[j];
        acc += diff_masspositions * masses[i] /
               std::pow(diff_masspositions.norm(), 3);
      }
    }
    forces[j] = acc * masses[j] / (4 * M_PI);
  }
// **********************************************************************
  return forces;
}

/* SAM_LISTING_END_2 */

/* SAM_LISTING_BEGIN_3 */
StarQuadTree::StarQuadTree(const std::vector<Eigen::Vector2d> &starpos,
                           const std::vector<double> &starmasses)
    : n(starpos.size()), starpos_(starpos), starmasses_(starmasses) {
  std::vector<unsigned int> star_idx(n);
  assertm(starmasses_.size() == n, "Size mismatch of star masses sequence");
  for (int i = 0; i < n; ++i) {
    star_idx[i] = i;
    assertm(((starpos_[i][0] >= 0.0) and (starpos_[i][0] <= 1.0) and
             (starpos_[i][1] >= 0.0) and (starpos_[i][1] <= 1.0)),
            "stars must be located insided unit square");
  }

  // Unit square is the bounding box for the whole group of stars
  Eigen::Matrix2d mainBbox;
  mainBbox << 0.0, 1.0, 0.0, 1.0;  // First row: x-coords, second row: y-coords

  root_ = std::move(std::make_unique<StarQuadTree::StarQuadTreeNode>(
      star_idx, mainBbox, *this));
}
/* SAM_LISTING_END_3 */

/* SAM_LISTING_BEGIN_X */
StarQuadTree::StarQuadTreeNode::StarQuadTreeNode(
    std::vector<unsigned int> star_idx, Eigen::Matrix2d bbox,
    StarQuadTree &tree)
    : star_idx_(star_idx), bbox_(bbox), mass(0.0), center({0, 0}) {
  assertm(star_idx_.size() > 0, "Can't create a node without stars...");
// **********************************************************************
// Code to be supplemented
for (unsigned int i{0}; i < star_idx_.size(); i++) {
  const Eigen::Vector2d starpos (tree.starpos_[star_idx_[i]]); // Initializes the pos for each star in each cluster
  mass += tree.starmasses_[star_idx_[i]];
  center += tree.starmasses_[star_idx_[i]] * starpos;
}
// Center of gravity for cluster of star associated with current node
center = center /mass;

if (star_idx_.size() == 1) {
  tree.no_leaves_++;
  return;
}
tree.no_clusters_++;

const double m1 = 0.5 * (bbox_(0, 0) + bbox_(0, 1));
const double m2 = 0.5 * (bbox_(1, 0) + bbox_(1, 1));

std::array<Eigen::Matrix2d, 4> sub_boxes;
sub_boxes[0] << bbox_(0, 0), m1, bbox(1, 0), m2;
sub_boxes[1] << bbox_(0, 0), m1, m2, bbox_(1, 1);
sub_boxes[2] << m1, bbox_(0, 1), bbox_(1, 0), m2;
sub_boxes[3] << m1, bbox_(0, 1), m2, bbox_(1, 1);
// Indes sets for son clusters:
std::array<std::vector<unsigned int>, 4> sub_indices;

// Check, which son box the stars lie in, cf.
for (unsigned int i {0}; i < star_idx_.size(); ++i) {
  if (tree.starpos_[star_idx_[i]][0] < m1) {   // left part
    if (tree.starpos_[star_idx_[i]][1] < m2) { // left bottom
      sub_indices[0].push_back(star_idx_[i]);
    } else {
      sub_indices[1].push_back(star_idx_[i]);  // left top
    }
  } else {                                     // right part
    if (tree.starpos_[star_idx_[i]][1] < m2) {
      sub_indices[2].push_back(star_idx_[i]);  // right bottom
    } else {
      sub_indices[3].push_back(star_idx_[i]);  // right top
    }
  }
}

// Recursion: son nodes are created if non-empty!
for (int i : {0, 1, 2, 3}) {
  if (!sub_indices[i].empty()) {
    sons_[i] = std::move(std::make_unique<StarQuadTree::StarQuadTreeNode>(sub_indices[i], sub_boxes[i], tree));
  }
}
// **********************************************************************
}

/* SAM_LISTING_END_X */

void StarQuadTree::outputQuadTree(std::ostream &o) const {
  o << "StarQuadTree with " << n << " stars, " << no_clusters_ << " clusters, "
    << no_leaves_ << " leaves" << std::endl;
  std::function<void(const StarQuadTreeNode *, int)> recursive_output =
      [&](const StarQuadTreeNode *node, int level) -> void {
    if (node != nullptr) {  // Safeguard
      const std::string indent(level, '>');
      o << indent << " Node on level " << level << ", star_idx_ = ";
      for (unsigned int idx : node->star_idx_) {
        o << idx << ' ';
      }
      o << std::endl;
      // indent << " bbox = " << node->bbox_ << std::endl
      o << indent << " center = " << node->center.transpose()
        << ", mass = " << node->mass << std::endl;
      for (unsigned int i : {0, 1, 2, 3}) {
        if (node->sons_[i]) {
          o << indent << " ### SON " << i << " ####" << std::endl;
          recursive_output((node->sons_[i]).get(), level + 1);
        }
      }
    }
  };
  recursive_output(root_.get(), 0);
}

StarQuadTreeClustering::StarQuadTreeClustering(
    const std::vector<Eigen::Vector2d> &starpos,
    const std::vector<double> &starmasses)
    : StarQuadTree(starpos, starmasses) {}

/* SAM_LISTING_BEGIN_5 */
bool StarQuadTreeClustering::isAdmissible(const StarQuadTreeNode &node,
                                          Eigen::Vector2d p, double eta) const {
  // Implements admissibility condition \lref{eq:admstar}
  bool admissible;
// **********************************************************************
// Code to be supplemented for 2-1d:
// Diameter:
  const double diam = (node.bbox_.col(0) - node.bbox_.col(1)).norm();

  auto intvdist = [](double a, double b, double x) -> double {
    if (b < a) std::swap(a, b);
    if (x < a) return (a - x);
    if (x > b) return (x - b);
    return 0.0;
  };
  const double dx = intvdist(node.bbox_(0, 0), node.bbox_(0, 1), p[0]);
  const double dy = intvdist(node.bbox_(1, 0), node.bbox_(1, 1), p[1]);
  const double dist = Eigen::Vector2d(dx, dy).norm();

  admissible = (dist > eta * diam);
// **********************************************************************
  return admissible;
}
/* SAM_LISTING_END_5 */

/* SAM_LISTING_BEGIN_6 */
Eigen::Vector2d StarQuadTreeClustering::forceOnStar(unsigned int j,
                                                    double eta) const {
  Eigen::Vector2d acc;  // For summation of force
  acc.setZero();
// **********************************************************************
// Code to be supplemented for 2-1e:
// **********************************************************************
  return acc;
}
/* SAM_LISTING_END_6 */

/* SAM_LISTING_BEGIN_7 */
std::vector<double> forceError(const StarQuadTreeClustering &qt,
                               const std::vector<double> &etas) {
  std::vector<double> error(etas.size());  // For returning errors
// **********************************************************************
// Code to be supplemented for 2-1f:
// **********************************************************************
  return error;
}
/* SAM_LISTING_END_7 */

/* SAM_LISTING_BEGIN_8 */
std::pair<double, double> measureRuntimes(unsigned int n, unsigned int n_runs) {
  assertm((n > 1), "At least two stars required!");
  // Initialize star positions
  std::vector<Eigen::Vector2d> pos = GravitationalForces::initStarPositions(n);
  // All stars have equal (unit) mass
  std::vector<double> mass(n, 1.0);
  double ms_exact = std::numeric_limits<double>::max();    // Time measured for exact evaluation
  double ms_cluster = std::numeric_limits<double>::max();  // Time taken for clustering-based evaluatiion
// **********************************************************************
// Code to be supplemented
// **********************************************************************
  return {ms_exact, ms_cluster};
}
/* SAM_LISTING_END_8 */

}  // namespace GravitationalForces

/* SAM_LISTING_BEGIN_9 */
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
/* SAM_LISTING_END_9 */
