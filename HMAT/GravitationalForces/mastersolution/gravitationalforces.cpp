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
          "Mismatch of sizes of masspositions and masses vectors");

  std::vector<Eigen::Vector2d> forces(n);

  Eigen::Vector2d acc;
  Eigen::Vector2d diff_masspositions;
  for (unsigned int j = 0; j < n; j++) {
    // Compute force on $\cob{j}$-th star
    acc.setZero();
    for (unsigned int i = 0; i < n; i++) {
      if (i != j) {
        diff_masspositions = masspositions[i] - masspositions[j];
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

  for (int i = 0; i < n; ++i) {
    star_idx[i] = i;
  }

  // @ToDo main box is hard-coded for now. A struct could be made for a box.
  Eigen::Matrix2d mainBbox;
  mainBbox << 0.0, 1.0, 0.0, 1.0;  // First row: x-axis, Second row: y-axis

  root_ = std::move(std::make_unique<StarQuadTree::StarQuadTreeNode>(
      star_idx, mainBbox, *this));
}
/* SAM_LISTING_END_3 */

/* SAM_LISTING_BEGIN_4 */
StarQuadTree::StarQuadTreeNode::StarQuadTreeNode(
    std::vector<unsigned int> star_idx, Eigen::Matrix2d bbox,
    const StarQuadTree &tree)
    : star_idx_(star_idx), bbox_(bbox), mass(0.0), center({0, 0}) {
  assertm(star_idx_.size() > 0, "Can't create a node without stars...");

  // Total mass
  for (unsigned int i = 0; i < star_idx_.size(); i++) {
    mass += tree.starmasses_[star_idx_[i]];
    center += tree.starmasses_[star_idx_[i]] * tree.starpos_[star_idx_[i]];
  }

  // Leaf
  if (star_idx_.size() == 1) return;

  // Center of gravity for non-leaf nodes
  center = center / mass;

  // Sub-boxes
  const double m1 = 0.5 * (bbox_(0, 0) + bbox_(0, 1));
  const double m2 = 0.5 * (bbox_(1, 0) + bbox_(1, 1));

  std::array<Eigen::Matrix2d, 4> sub_boxes;

  sub_boxes[0] << bbox_(0, 0), m1, bbox_(1, 0), m2;
  sub_boxes[1] << bbox_(0, 0), m1, m2, bbox_(1, 1);
  sub_boxes[2] << m1, bbox_(0, 1), bbox_(1, 0), m2;
  sub_boxes[3] << m1, bbox_(0, 1), m2, bbox_(1, 1);

  // Sub-indices
  std::array<std::vector<unsigned int>, 4> sub_indices;

  //@ToDo Bother with sorting?
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

  // Sons get initialized if necessary, otherwise they stay as nullptr
  for (int i : {0, 1, 2, 3}) {
    if (!sub_indices[i].empty())
      sons_[i] = std::move(std::make_unique<StarQuadTree::StarQuadTreeNode>(
          sub_indices[i], sub_boxes[i], tree));
  }
}
/* SAM_LISTING_END_4 */

StarQuadTreeClustering::StarQuadTreeClustering(
    const std::vector<Eigen::Vector2d> &starpos,
    const std::vector<double> &starmasses)
    : StarQuadTree(starpos, starmasses) {}

/* SAM_LISTING_BEGIN_5 */
// TODO: Implementation of distance is not correct
bool StarQuadTreeClustering::isAdmissible(const StarQuadTreeNode &node,
                                          Eigen::Vector2d p, double eta) const {
  if (node.isLeaf()) return true;  // Leafs are admissible

  // @ToDo This line is called four times, once for each brother with same diam.
  // Make a struct for bbox to avoid multiple computations of same diam.
  double diam = (node.bbox_.col(1) - node.bbox_.col(0)).norm();

  // Supremum distance between p and the four corner points of bbox
  double dist = (p - node.bbox_.col(0)).norm();  // Corner point 1
  double tmp = (p - node.bbox_.col(1)).norm();
  if (tmp > dist) dist = tmp;  // Corner point 2
  tmp = (p - node.bbox_.diagonal()).norm();
  if (tmp > dist) dist = tmp;  // Corner point 3
  tmp = (p - Eigen::Vector2d{node.bbox_(0, 1), node.bbox_(1, 0)}).norm();
  if (tmp > dist) dist = tmp;  // Corner point 4

  return (dist > eta * diam);  // Admissibility condition

}
/* SAM_LISTING_END_5 */

/* SAM_LISTING_BEGIN_6 */
Eigen::Vector2d StarQuadTreeClustering::forceOnStar(unsigned int j,
                                                    double eta) const {
  Eigen::Vector2d acc;
  acc.setZero();
  Eigen::Vector2d diff_masspositions;

  // Trick: Recursive lambda function !
  std::function<void(const std::unique_ptr<StarQuadTreeNode> &)> traverse =
      [&](const std::unique_ptr<StarQuadTreeNode> &node) {
        if (!node) return;  // In case if the cluster has no stars in it
        // @ToDo check this with boss. Point p is considered not admissible, if
        // it lies itself in the current subcluster.
        // Alternatively, check at least the following condition:
        // node->star_idx_[0] != j (if node is a Leaf, its star cannot be x^j!)
        if (!std::count(node->star_idx_.begin(), node->star_idx_.end(), j) &&
            isAdmissible(*node, starpos_[j], eta)) {
          diff_masspositions = node->center - starpos_[j];
          acc += diff_masspositions * node->mass /
                 pow(diff_masspositions.norm(), 3);
        } else {  // traverse further, if current sub-cluster is not admissible
          traverse(node->sons_[0]);
          traverse(node->sons_[1]);
          traverse(node->sons_[2]);
          traverse(node->sons_[3]);
        }
      };
  // Start of recursion
  traverse(this->root_);
  // Multiply with forefactor in \prbeqref{eq:fjs}
  return (acc * starmasses_[j] / (4 * M_PI));
}
/* SAM_LISTING_END_6 */

/* SAM_LISTING_BEGIN_7 */
std::vector<double> forceError(const StarQuadTreeClustering &qt,
                               const std::vector<double> &etas) {
  std::vector<double> error(etas.size());
  std::vector<Eigen::Vector2d> exact_forces;

  // time the evaluation of computeForces_direct
  if constexpr (TIMING_MODE > 0) {
    auto t1 = std::chrono::high_resolution_clock::now();
    for (unsigned int trange = 0; trange < TIMING_MODE; trange++) {
      exact_forces = computeForces_direct(qt.starpos_, qt.starmasses_);
    }
    auto t2 = std::chrono::high_resolution_clock::now();
    /* Getting number of milliseconds as a double. */
    std::chrono::duration<double, std::milli> ms_double =
        (t2 - t1) / TIMING_MODE;
    std::cout << "Runtime computeForces_direct= " << ms_double.count()
              << "ms\n";
  } else {  // evaluate computeForces_direct without timings
    exact_forces = computeForces_direct(qt.starpos_, qt.starmasses_);
  }

  std::vector<double> normed_errors(qt.n);
  for (int eta_i = 0; eta_i < etas.size(); eta_i++) {
    if constexpr (TIMING_MODE > 0) {  // compute and time forceOnStar
      auto t1 = std::chrono::high_resolution_clock::now();
      for (unsigned int trange = 0; trange < TIMING_MODE; trange++) {
        for (unsigned int j = 0; j < qt.n; j++) {
          normed_errors[j] =
              (exact_forces[j] - qt.forceOnStar(j, etas[eta_i])).norm();
        }
      }
      auto t2 = std::chrono::high_resolution_clock::now();
      /* Getting number of milliseconds as a double. */
      std::chrono::duration<double, std::milli> ms_double =
          (t2 - t1) / TIMING_MODE;
      std::cout << "Runtime forceOnStar[eta=" << etas[eta_i]
                << "]= " << ms_double.count() << "ms\n";
    } else {  // compute forceOnStar without timings
      for (unsigned int j = 0; j < qt.n; j++) {
        normed_errors[j] =
            (exact_forces[j] - qt.forceOnStar(j, etas[eta_i])).norm();
        // std::cout<<"Normed error "<<j<<": "<<normed_errors[j]<<std::endl;
      }
    }

    error[eta_i] =
        *std::max_element(normed_errors.begin(), normed_errors.end());
  }

  return error;
}
/* SAM_LISTING_END_7 */

/* SAM_LISTING_BEGIN_Y */
std::pair<double, double> measureRuntimes(unsigned int n) {
  // TO BE SUPPLEMENTED
}
/* SAM_LISTING_END_Y */

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
