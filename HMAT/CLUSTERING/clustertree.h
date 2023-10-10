/***********************************************************************
 *                                                                     *
 * Code for Course "Advanced Numerical Methods for CSE"                *
 * (Prof. Dr. R. Hiptmair)                                             *
 * Author: R.H.                                                        *
 * Date: Nov 18, 2017                                                  *
 * (C) Seminar for Applied Mathematics, ETH Zurich                     *
 * This code can be freely used for non-commercial purposes as long    *
 * as this header is left intact.                                      *
 ***********************************************************************/
#ifndef CLUSTERTREE_H_
#define CLUSTERTREE_H_

// General includes
#include <Eigen/Dense>
#include <algorithm>
#include <array>
#include <cmath>
#include <exception>
#include <iostream>
#include <limits>
#include <utility>
#include <vector>

namespace HMAT {
/** @brief Data structure for a collocation point
 A collocation point has an index and coordinates
 @tparam DIM dimension of ambient space
*/
/* SAM_LISTING_BEGIN_1 */
template <int DIM>  // dimension \cob{$d$} as template argument
struct Point {
  constexpr static std::size_t dim = DIM;
  std::size_t idx;                  // number of collocation point
  Eigen::Matrix<double, DIM, 1> x;  // coordinate vector
};
/* SAM_LISTING_END_1 */

/** @brief Data structure for bounding box
 *  a bounding box is defined by two position vectors corresponding to
 * lower-left and upper-right corners, defining d intervals
 *  @tparam d dimension of ambient spacew
 */
/* SAM_LISTING_BEGIN_2 */
template <int DIM>  // dimension \cob{$d$} as template argument
struct BBox {
  constexpr static std::size_t dim = DIM;
  // Bounding box from sequence of points
  explicit BBox(const std::vector<Point<DIM>> pts);
  // Size \cob{$\diam(B)$} of a bounding box
  [[nodiscard]] double diam() const {
    return (maxc - minc).cwiseAbs().maxCoeff();
  }
  // Coordinate vectors of Corner points of bounding box
  Eigen::Matrix<double, DIM, 1> minc;  // Lower-left corner
  Eigen::Matrix<double, DIM, 1> maxc;  // Upper-right corner
};
/* SAM_LISTING_END_2 */

// distance of 1D intervals \cob{$\cintv{a,b}$} and \cob{$\cintv{c,d}$}
double dist(double a, double b, double c, double d);

// distance of d-dimensional boxes
template <int DIM>
double dist(const BBox<DIM> &bx, const BBox<DIM> &by) {
  double dst = 0.0;
  for (int l = 0; l < DIM; ++l) {
    dst += pow(dist(bx.minc[l], bx.maxc[l], by.minc[l], by.maxc[l]), 2);
  }
  return sqrt(dst);
}

// output bounding box
template <int DIM>
std::ostream &operator<<(std::ostream &o, const BBox<DIM> &box) {
  return o << "BBOX: " << box.minc.transpose() << " - " << box.maxc.transpose()
           << ' ';
}

/* SAM_LISTING_BEGIN_3 */
template <int DIM>
BBox<DIM>::BBox(const std::vector<Point<DIM>> pts) {
  double tmp;
  minc = std::numeric_limits<double>::max() *
         Eigen::Matrix<double, DIM, 1>::Ones();
  maxc = -minc;
  for (const Point<DIM> &v : pts) {
    for (int l = 0; l < DIM; l++) {
      const double tmp = v.x[l];
      if (tmp < minc[l]) {
        minc[l] = tmp;
      }
      if (tmp > maxc[l]) {
        maxc[l] = tmp;
      }
    }
  }
}  // end BBox constructor
/* SAM_LISTING_END_3 */

/** @brief Data structure for the node of a binary cluster tree
   A node of a cluster tree contains a set of collocation points.
   @tparam DIM dimension of ambient space
*/
/* SAM_LISTING_BEGIN_4 */
template <int DIM>
class CtNode {
 public:
  constexpr static std::size_t dim = DIM;
  // Constructors taking a sequence of points
  explicit CtNode(const std::vector<Point<DIM>> _pts, int _dir = 0)
    : pts(std::move(_pts)), sons{nullptr, nullptr}, dir(_dir),
      clust_sect_vec(_pts.size()) {}
  // Destructor (also attempts to destroy sons!)
  virtual ~CtNode() {
    if (sons[0]) {
      delete sons[0];
    }
    if (sons[1]) {
      delete sons[1];
    }
  }
  // Number of indices owned by the cluster
  [[nodiscard]] std::size_t noIdx() const { return pts.size(); }
  // Function \cob{$\Ci$}: access to owned indices
  [[nodiscard]] std::vector<std::size_t> I() const;
  // Access to bounding box (computed on the fly)
  [[nodiscard]] BBox<DIM> getBBox() const { return BBox<DIM>(pts); }
  // Is the node a leaf node ?
  [[nodiscard]] virtual bool isLeaf() const {
    return (!(sons[0]) and !(sons[1]));
  }
  // Output operator or recursive output
  template <int dim>
  friend std::ostream &operator<<(std::ostream &o, const CtNode<dim> &node);
  // Public data member: Pointers to two (binary tree!) sons
  std::array<CtNode *, 2> sons;
  // Public data member: Points contained in the cluster
  std::vector<Point<DIM>> pts;
  // Temporary storage for cluster-associated vector sections
  Eigen::VectorXd clust_sect_vec;  
  // Direction for sorting, passed by the constructor
  int dir;
};
/* SAM_LISTING_END_4 */

/* SAM_LISTING_BEGIN_8 */
template <int DIM>
std::vector<size_t> CtNode<DIM>::I() const {
  std::vector<size_t> idx;
  for (const Point<DIM> &pt : pts) idx.push_back(pt.idx);
  return idx;
}
/* SAM_LISTING_END_8 */

  // Output operator for a plain cluster
template <int dim>
std::ostream &operator<<(std::ostream &o, const CtNode<dim> &node) {
  o << "Node: " << node.pts.size() << " points, ";
  if (node.isLeaf()) {
    o << "LEAF ";
  }
  o << node.getBBox() << ": ";
  for (const Point<dim> &pt : node.pts) {
    o << "[";
    for (int l = 0; l < dim; l++) {
      o << pt.x[l] << ' ';
    }
    o << "]";
  }
  o << std::endl;
  if (node.sons[0]) {
    o << *(node.sons[0]);
  }
  if (node.sons[1]) {
    o << *(node.sons[1]);
  }
  return o;
}

// For debugging purposes: non-recursive output for 1D
void outCTNode(const CtNode<1> &ipnode, std::ostream &o = std::cout);

  
/** @brief Data structure for a cluster tree.
    A cluster tree builds and describes a multilevel partitioning of a set of
   points
   @tparam Node Data type for a single node, see CtNode
*/
/* SAM_LISTING_BEGIN_5 */
template <class NODE>
class ClusterTree {
 public:
  using node_t = NODE;
  constexpr static std::size_t dim = NODE::dim;  // space dimension d
  // Idle constructor
  ClusterTree() : root(nullptr) {}
  // Effective Constructor taking a sequence of points
  // (needed, because polynorphism not supported in constructor)
  void init(const std::vector<Point<dim>> &pts, std::size_t minpts = 1);
  // Recursive destruction
  virtual ~ClusterTree() {
    if (root) delete root;
  }
  // Output of tree
  template <class Nd>
  friend std::ostream &operator<<(std::ostream &o, const ClusterTree<Nd> &T);

 protected:
  // Recursive construction
  virtual void buildRec(NODE *nptr, std::size_t minpts);
  // Node factory
  virtual NODE *createNode(const std::vector<Point<dim>> pts, int dir) {
    return new NODE(pts, dir);
  }

 public:
  NODE *root;  // pointer to root node
};
/* SAM_LISTING_END_5 */

/* SAM_LISTING_BEGIN_6 */
template <class NODE>
void ClusterTree<NODE>::init(const std::vector<Point<dim>> &pts,
                             std::size_t minpts) {
  if (!(root = createNode(pts, 0))) {
    throw(std::runtime_error("Cannot allocate root"));
  }
  if (minpts < 1) {
    throw(std::runtime_error("minpts must be at least 1"));
  }
  buildRec(root, minpts);
}
/* SAM_LISTING_END_6 */

template <class NODE>
std::ostream &operator<<(std::ostream &o, const ClusterTree<NODE> &T) {
  return o << "ROOT: " << *(T.root) << "END CLUSTER TREE" << std::endl;
}

/* SAM_LISTING_BEGIN_7 */
template <class NODE>
void ClusterTree<NODE>::buildRec(NODE *nptr, std::size_t minpts) {
  const std::size_t n = nptr->pts.size();  // Number of held indices
  // Leaf, if minimal number of indices reached
  if (n > minpts) {  // \Label[line]{brc:1}
    // Points have to be copied and sorted according to direction dir
    std::vector<Point<dim>> tpts(nptr->pts);
    // next sorting direction
    const int dir = (nptr->dir + 1) % dim;
    // call sort function from standard library
    std::sort(tpts.begin(), tpts.end(),
              [dir](const Point<dim> &p1, const Point<dim> &p2) -> bool {
                return (bool)(p1.x[dir] < p2.x[dir]);
              });
    // Split point sequence and construct sons
    const std::size_t m = n / 2;  // integer arithmeric, m>0 ensured
    const std::vector<Point<dim>> low_pts(tpts.cbegin(), tpts.cbegin() + m);
    // First son gets ``lower half'' of sorted points
    if (!(nptr->sons[0] = createNode(low_pts, dir))) {
      throw(std::runtime_error("Cannot allocate first son"));
    }
    buildRec(nptr->sons[0], minpts);  // recurse into first son
    // Second son get ``upper half'' of sorted points
    const std::vector<Point<dim>> up_pts(tpts.cbegin() + m, tpts.cend());
    if (!(nptr->sons[1] = createNode(up_pts, dir))) {
      throw(std::runtime_error("Cannot allocate second son"));
    }
    buildRec(nptr->sons[1], minpts);  // recurse into 2nd son
  }
}
/* SAM_LISTING_END_7 */

}  // end namespace HMAT

#endif
