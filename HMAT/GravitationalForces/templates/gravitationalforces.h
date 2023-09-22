/**
 * @file gravitationalforces.h
 * @brief NPDE homework GravitationalForces code
 * @author
 * @date
 * @copyright Developed at SAM, ETH Zurich
 */

#ifndef GF_H_
#define GF_H_

#include <Eigen/Dense>
#include <cassert>
#include <iomanip>
#include <iostream>
#include <memory>
#include <ostream>
#include <utility>
#include <vector>
#include <array>
#define assertm(exp, msg) assert(((void)msg, exp))

namespace GravitationalForces {

/** @brief Random distribution of mass points in the unit square;
 *
 * @param n number of masspoints
 * @param mindist minimal distance of two mass points
 *
 * A minimal distance is imposed, smaller than \f$ \frac{1}{2\sqrt{2}}\f$
 */
std::vector<Eigen::Vector2d> initStarPositions(unsigned int n,
                                               double mindist = 1.0);

/** @brief Exact computation of forces on all stars
 *
 * @param masspositions sequence of coordinate vectors of stars treated as mass
 * points
 * @param masses sequence of masses of stars
 */
std::vector<Eigen::Vector2d> computeForces_direct(
    const std::vector<Eigen::Vector2d> &masspositions,
    const std::vector<double> &masses);

/** @brief Class representing a quadtree of "equivalent stars"
 *
 */
/* SAM_LISTING_BEGIN_1 */
class StarQuadTree {
 public:
  StarQuadTree() = delete;
  StarQuadTree(const StarQuadTree &) = delete;
  StarQuadTree(StarQuadTree &&) noexcept = default;
  StarQuadTree &operator=(const StarQuadTree &) = delete;
  StarQuadTree &operator=(StarQuadTree &&) noexcept = default;
  // Constructor, which also builds the tree
  StarQuadTree(const std::vector<Eigen::Vector2d> &starpos,
               const std::vector<double> &starmasses);
  virtual ~StarQuadTree() = default;

 protected:
  struct StarQuadTreeNode {
    // Constructor: does the recursive construction of the quadtree
    StarQuadTreeNode(std::vector<unsigned int> star_idx, Eigen::Matrix2d bbox,
                     StarQuadTree &tree);
    virtual ~StarQuadTreeNode() = default;
    // In this code: leaf nodes hold only a single star!
    [[nodiscard]] inline bool isLeaf() const {
      return this->star_idx_.size() == 1;
    }
    std::vector<unsigned int> star_idx_;  // Indices of stars in sub-cluster
    Eigen::Matrix2d bbox_;                // Bounding box of sub-cluster
    Eigen::Vector2d center;               // Center of gravity of sub-cluster
    double mass;                          // Total mass of stars in sub-cluster
    // Pointers to children nodes
    std::array<std::unique_ptr<StarQuadTreeNode>, 4> sons_{nullptr};
  };

 public:
  const int n;                                  // Total number of stars
  std::unique_ptr<StarQuadTreeNode> root_;      // Root node
  const std::vector<Eigen::Vector2d> starpos_;  // "Point star" positions
  const std::vector<double> starmasses_;        // Star masses
  unsigned int no_leaves_{0};                   // number of leaves
  unsigned int no_clusters_{0};  // number of genuine star clusters

  // Printing the tree (useful only for small trees)
  void outputQuadTree(std::ostream &o) const;
};
/* SAM_LISTING_END_1 */

/* SAM_LISTING_BEGIN_2 */
class StarQuadTreeClustering : public StarQuadTree {
 public:
  StarQuadTreeClustering() = delete;
  StarQuadTreeClustering(const StarQuadTreeClustering &) = delete;
  StarQuadTreeClustering(StarQuadTreeClustering &&) noexcept = default;
  StarQuadTreeClustering &operator=(const StarQuadTreeClustering &) = delete;
  StarQuadTreeClustering &operator=(StarQuadTreeClustering &&) noexcept =
      default;
  // The only non-trivial constructor
  StarQuadTreeClustering(const std::vector<Eigen::Vector2d> &starpos,
                         const std::vector<double> &starmasses);
  virtual ~StarQuadTreeClustering() = default;

  // Decide the admissibility of a sub-cluster with respect to a point
  [[nodiscard]] bool isAdmissible(const StarQuadTreeNode &node,
                                  Eigen::Vector2d p, double eta) const;
  // Compute total graviational force on a single star with index j
  [[nodiscard]] Eigen::Vector2d forceOnStar(unsigned int j, double eta) const;
};
/* SAM_LISTING_END_2 */

std::vector<double> forceError(const StarQuadTreeClustering &qt,
                               const std::vector<double> &etas);

// Runtime measurements
std::pair<double, double> measureRuntimes(unsigned int n,
                                          unsigned int n_runs = 5);

}  // namespace GravitationalForces

// Code demonstrating how to measure runtimes
void runtimemeasuredemo(void);

#endif
