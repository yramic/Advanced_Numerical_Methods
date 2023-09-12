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
          "Mismatch of sizes of masspositions and masses vectors");

  std::vector<Eigen::Vector2d> forces(n);

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
            "stars must be located insuided unit square");
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

  // Total mass
  for (unsigned int i = 0; i < star_idx_.size(); i++) {
    const Eigen::Vector2d starpos(tree.starpos_[star_idx_[i]]);
    mass += tree.starmasses_[star_idx_[i]];
    center += tree.starmasses_[star_idx_[i]] * starpos;
    assertm(((starpos[0] >= bbox(0, 0)) and (starpos[0] <= bbox(0, 1)) and
             (starpos[1] >= bbox(1, 0)) and (starpos[1] <= bbox(1, 1))),
            "Star out of bounding box!");
  }
  // Center of gravity for cluster of star associated with current node
  center = center / mass;

  // In this implementation leaf nodes contain a single star, whose position
  // also agrees with the center of the node cluster.
  if (star_idx_.size() == 1) {
    tree.no_leaves_++;
    return;
  }
  tree.no_clusters_++;
  // Sub-boxes obtained by halfing bounding box in each direction
  const double m1 = 0.5 * (bbox_(0, 0) + bbox_(0, 1));
  const double m2 = 0.5 * (bbox_(1, 0) + bbox_(1, 1));
  std::array<Eigen::Matrix2d, 4> sub_boxes;
  sub_boxes[0] << bbox_(0, 0), m1, bbox_(1, 0), m2;
  sub_boxes[1] << bbox_(0, 0), m1, m2, bbox_(1, 1);
  sub_boxes[2] << m1, bbox_(0, 1), bbox_(1, 0), m2;
  sub_boxes[3] << m1, bbox_(0, 1), m2, bbox_(1, 1);
  // Index sets for son clusters
  std::array<std::vector<unsigned int>, 4> sub_indices;

  // Check, which son box the stars lie in, cf. \lref{pc:geoclust}
  for (unsigned int i = 0; i < star_idx_.size(); i++) {
    if (tree.starpos_[star_idx_[i]][0] < m1) {    // left part
      if (tree.starpos_[star_idx_[i]][1] < m2) {  // left bottom
        sub_indices[0].push_back(star_idx_[i]);
      } else {  // left top
        sub_indices[1].push_back(star_idx_[i]);
      }
    } else {                                      // right part
      if (tree.starpos_[star_idx_[i]][1] < m2) {  // right bottom
        sub_indices[2].push_back(star_idx_[i]);
      } else {  // right top
        sub_indices[3].push_back(star_idx_[i]);
      }
    }
  }
  // Recursion: son nodes are created, if non-empty
  for (int i : {0, 1, 2, 3}) {
    if (!sub_indices[i].empty())
      sons_[i] = std::move(std::make_unique<StarQuadTree::StarQuadTreeNode>(
          sub_indices[i], sub_boxes[i], tree));
  }
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
  // Diameter of bounding box is the distance of its opposite corners
  const double diam = (node.bbox_.col(0) - node.bbox_.col(1)).norm();
  // To compute the distance of a point from an axis-aligned box we have to
  // distinguish nine different cases, which can conveniently be done by first
  // introducing a function computing distances in 1D.
  auto intvdist = [](double a, double b, double x) -> double {
    if (b < a) std::swap(a, b);
    if (x < a) return (a - x);
    if (x > b) return (x - b);
    return 0.0;
  };
  const double dx = intvdist(node.bbox_(0, 0), node.bbox_(0, 1), p[0]);
  const double dy = intvdist(node.bbox_(1, 0), node.bbox_(1, 1), p[1]);
  const double dist = Eigen::Vector2d(dx, dy).norm();
  return (dist > eta * diam);  // Admissibility condition \lref{eq:admstar}
}
/* SAM_LISTING_END_5 */

/* SAM_LISTING_BEGIN_6 */
Eigen::Vector2d StarQuadTreeClustering::forceOnStar(unsigned int j,
                                                    double eta) const {
  Eigen::Vector2d acc;  // For summation of force
  acc.setZero();

  // Trick: Recursive lambda function capturing the whole object
  std::function<void(const StarQuadTreeNode *)> traverse =
      [&](const StarQuadTreeNode *node) {
        if (node == nullptr)
          return;  // In case if the cluster has no stars in it
        // Since leaf clusters contain only a single star we can treat them in
        // the same way as admissible clusters.
        if (isAdmissible(*node, starpos_[j], eta) or (node->isLeaf())) {
          const Eigen::Vector2d diff_masspositions = node->center - starpos_[j];
          const double dist = diff_masspositions.norm();
          if (dist > 1.0E-10) {
            acc += diff_masspositions * node->mass / pow(dist, 3);
          }
        } else {  // traverse further, if current sub-cluster is not
                  // admissible
          traverse(node->sons_[0].get());
          traverse(node->sons_[1].get());
          traverse(node->sons_[2].get());
          traverse(node->sons_[3].get());
        }
      };
  // Start of recursion
  traverse(this->root_.get());
  // Multiply with forefactor in \prbeqref{eq:fjs}
  return (acc * starmasses_[j] / (4 * M_PI));
}
/* SAM_LISTING_END_6 */

/* SAM_LISTING_BEGIN_7 */
std::vector<double> forceError(const StarQuadTreeClustering &qt,
                               const std::vector<double> &etas) {
  std::vector<double> error(etas.size());  // For returning errors
  std::vector<Eigen::Vector2d> exact_forces{
      computeForces_direct(qt.starpos_, qt.starmasses_)};

  std::vector<double> normed_errors(qt.n);
  // Compute errors for different values of the admissibility parameter
  for (int eta_i = 0; eta_i < etas.size(); eta_i++) {
    for (unsigned int j = 0; j < qt.n; j++) {
      const Eigen::Vector2d force_j = qt.forceOnStar(j, etas[eta_i]);
      normed_errors[j] = (exact_forces[j] - force_j).norm();
    }
    error[eta_i] =
        *std::max_element(normed_errors.begin(), normed_errors.end());
  }
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
  // Build quad tree of stars
  StarQuadTreeClustering qt(pos, mass);
  // Admissibility parameter
  const double eta = 1.5;

  // Runtime for exact computation of forces with effort $O(n^2)$
  std::vector<Eigen::Vector2d> forces(n);
  double ms_exact = 0.0;
  for (int r = 0; r < n_runs; ++r) {
    auto t1_exact = std::chrono::high_resolution_clock::now();
    forces = computeForces_direct(qt.starpos_, qt.starmasses_);
    auto t2_exact = std::chrono::high_resolution_clock::now();
    /* Getting number of milliseconds as a double. */
    std::chrono::duration<double, std::milli> ms_double = (t2_exact - t1_exact);
    ms_exact = std::max(ms_exact, ms_double.count());
  }
  std::cout << "n = " << n << " : runtime computeForces_direct= " << ms_exact
            << "ms\n";

  // Runtime for cluster-based approximate evaluation, cost $O(n\log n)$
  double ms_cluster = 0.0;
  for (int r = 0; r < n_runs; ++r) {
    auto t1_cluster = std::chrono::high_resolution_clock::now();
    for (unsigned int j = 0; j < qt.n; j++) {
      forces[j] = qt.forceOnStar(j, eta);
    }
    auto t2_cluster = std::chrono::high_resolution_clock::now();
    /* Getting number of milliseconds as a double. */
    std::chrono::duration<double, std::milli> ms_double =
        (t2_cluster - t1_cluster);
    ms_cluster = std::max(ms_cluster, ms_double.count());
  }
  std::cout << "n = " << n << " : runtime forceOnStar[eta=" << eta
            << "]= " << ms_cluster << "ms\n";
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
