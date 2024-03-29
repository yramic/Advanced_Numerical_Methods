/***********************************************************************
 *                                                                     *
 * Code for Course "Advanced Numerical Methods for CSE"                *
 * (Prof. Dr. R. Hiptmair)                                             *
 * Author: R.H., Peiyuan Xie                                           *
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
#include <memory>
#include <stack>
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
   A node of a cluster tree contains index set, node number and offset.
   @tparam DIM dimension of ambient space
*/
/* SAM_LISTING_BEGIN_4 */
template <int DIM>
class CtNode {
 public:
  constexpr static std::size_t dim = DIM;
  // Constructors taking a sequence of indices, offset, node number and sorting
  // direaction
  explicit CtNode(const std::vector<size_t> _I, int _offset, int _nodeNumer,
                  int _dir = 0)
      : I(std::move(_I)),
        sons{nullptr, nullptr},
        dir(_dir),
        offset(_offset),
        nodeNumber(_nodeNumer) {}

  // Number of indices owned by the cluster
  [[nodiscard]] std::size_t noIdx() const { return I.size(); }
  // Is the node a leaf node ?
  [[nodiscard]] virtual bool isLeaf() const {
    return (!(sons[0]) and !(sons[1]));
  }
  // Public data member: Pointers to two (binary tree!) sons
  std::array<std::unique_ptr<CtNode>, 2>
      sons;  // smart pointer, no need for destructor
  // Public data member: Index set of the cluster
  std::vector<size_t> I;
  // Public data member: Offset indicating where local points are stored in
  // global vector
  int offset;
  // Public data member: Node number, to reference cluster-local data kept
  // elsewhere
  int nodeNumber;
  // Direction for sorting, passed by the constructor
  int dir;
};
/* SAM_LISTING_END_4 */

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
  ClusterTree(const std::vector<Point<dim>> &pts, std::size_t minpts = 1)
      : root(nullptr) {
    numNodes = 0;
    int offset = 0;
    // Append local collocation points to the end of global vector
    ptsT.insert(ptsT.end(), pts.begin(), pts.end());
    // Build index set for root
    std::vector<size_t> idx;
    for (const Point<dim> &pt : pts) idx.push_back(pt.idx);
    root = std::make_unique<NODE>(idx, offset, numNodes++, 0);
    offset += idx.size();
    if (!root) {
      throw(std::runtime_error("Cannot allocate root"));
    }
    if (minpts < 1) {
      throw(std::runtime_error("minpts must be at least 1"));
    }
    // Recursive construction of child nodes
    std::function<void(NODE * nptr)> buildrec = [&](NODE *nptr) -> void {
      const std::size_t n = nptr->noIdx();  // Number of held indices
      // Leaf, if minimal number of indices reached
      if (n > minpts) {  // \Label[line]{brc:1}
        // Points have to be copied and sorted according to direction dir
        std::vector<Point<dim>> tpts(
            ptsT.begin() + nptr->offset,
            ptsT.begin() + nptr->offset + nptr->noIdx());
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
        // Append local collocation points to the end of global vector
        ptsT.insert(ptsT.end(), low_pts.begin(), low_pts.end());
        // Build index set for child node
        idx.clear();
        for (const Point<dim> &pt : low_pts) idx.push_back(pt.idx);
        // First son gets ``lower half'' of sorted points
        nptr->sons[0] = std::make_unique<NODE>(idx, offset, numNodes++, dir);
        offset += idx.size();
        if (!nptr->sons[0]) {
          throw(std::runtime_error("Cannot allocate first son"));
        }
        buildrec(nptr->sons[0].get());  // recurse into first son
        // Append local collocation points to the end of global vector
        const std::vector<Point<dim>> up_pts(tpts.cbegin() + m, tpts.cend());
        ptsT.insert(ptsT.end(), up_pts.begin(), up_pts.end());
        // Build index set for child node
        idx.clear();
        for (const Point<dim> &pt : up_pts) {
          idx.push_back(pt.idx);
        }
        // Second son get ``upper half'' of sorted points
        nptr->sons[1] = std::make_unique<NODE>(idx, offset, numNodes++, dir);
        offset += idx.size();
        if (!nptr->sons[1]) {
          throw(std::runtime_error("Cannot allocate second son"));
        }
        buildrec(nptr->sons[1].get());  // recurse into 2nd son
      } 
    };
    buildrec(root.get());
  }

  // Access to bounding box of the given node (computed on the fly)
  virtual BBox<dim> getBBox(const NODE *nptr) const {
    return BBox<dim>(
        std::vector<Point<dim>>(ptsT.begin() + nptr->offset,
                                ptsT.begin() + nptr->offset + nptr->noIdx()));
  }

  // Output of tree
  template <class Nd>
  friend std::ostream &operator<<(std::ostream &o, const ClusterTree<Nd> &T);

 public:
  std::unique_ptr<NODE> root;    // pointer to root node
  std::vector<Point<dim>> ptsT;  // vector containing points for each node
  int numNodes;                  // total number of nodes in the tree
};
/* SAM_LISTING_END_5 */

// For debug, using iteration instead of recursion
template <class NODE>
std::ostream &operator<<(std::ostream &o, const ClusterTree<NODE> &T) {
  // o << "ROOT: ";
  o << "CLUSTER TREE" << std::endl;
  std::stack<NODE *> nodes;
  nodes.push(T.root.get());
  while (!nodes.empty()) {
    NODE *node = nodes.top();
    nodes.pop();
    o << "Node: " << node->noIdx() << " points, ";
    if (node->isLeaf()) {
      o << "LEAF ";
    }
    o << T.getBBox(node) << ": ";
    for (int i = 0; i < node->noIdx(); i++) {
      o << "[";
      for (int l = 0; l < NODE::dim; l++) {
        o << T.ptsT[node->offset + i].x[l] << ' ';
      }
      o << "]";
      // o << " " << node->nodeNumber << " " << node->offset;
    }
    o << std::endl;
    if (node->sons[1]) {
      nodes.push(node->sons[1].get());
    }
    if (node->sons[0]) {
      nodes.push(node->sons[0].get());
    }
  }
  // o << "END CLUSTER TREE" << std::endl;
  return o;
}

}  // end namespace HMAT

#endif
