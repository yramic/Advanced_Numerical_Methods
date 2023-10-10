/***********************************************************************
 *                                                                     *
 * Code for Course "Advanced Numerical Methods for CSE"                *
 * (Prof. Dr. R. Hiptmair)                                             *
 * Author: R.H.                                                        *
 * Date: August 2023
 * (C) Seminar for Applied Mathematics, ETH Zurich                     *
 * This code can be freely used for non-commercial purposes as long    *
 * as this header is left intact.                                      *
 ***********************************************************************/

#ifndef MATPART_H_
#define MATPART_H_

#include <cassert>
#include <cstddef>
#include <memory>

#include "clustertree.h"

#define assertm(exp, msg) assert(((void)msg, exp))

namespace HMAT {
/** @brief Data structure for both far-field and near-field blocks */
/* SAM_LISTING_BEGIN_9 */
template <class NODE>
struct IndexBlock {
  using node_t = NODE;
  constexpr static std::size_t dim = NODE::dim;
  // Constructors extracts indices from clusters
  IndexBlock(NODE &_nx, NODE &_ny)
      : nx(_nx), ny(_ny), i_idx(_nx.I), j_idx(ny.I) {}
  virtual ~IndexBlock() = default;
  NODE &nx, &ny;                           // contributing clusters
  const std::vector<size_t> i_idx, j_idx;  // contained indices
};
/* SAM_LISTING_END_9 */

/* SAM_LISTING_BEGIN_A */
template <class Tree>
class BlockPartition {
 public:
  using node_t = typename Tree::node_t;
  // Idle constructor
  BlockPartition(std::shared_ptr<Tree> _rowT, std::shared_ptr<Tree> _colT)
      : rowT(_rowT), colT(_colT) {
    assertm((rowT != nullptr), "No valid x-tree!");
    assertm((colT != nullptr), "No valid y-tree!");
  }
  // Trigger recursive construction of partition
  // (Needed, because polymorphic functions not available in constructor)
  void init(double eta0 = 0.5);
  virtual ~BlockPartition() = default;
  // Size of the matrix
  [[nodiscard]] size_t cols() const { return colT->root->noIdx(); }
  [[nodiscard]] size_t rows() const { return rowT->root->noIdx(); }
  // Admissibility condition \cob{$\adm$}, see \cref{def:ac}
  [[nodiscard]] virtual bool adm(const node_t *nx, const node_t *ny,
                                 double eta0) const;

 protected:
  // Recursive construction from cluster pair
  virtual void buildRec(node_t *nx, node_t *ny, double eta0);
  // Construct an instance of far-field block type
  virtual IndexBlock<node_t> makeFarFieldBlock(node_t &nx, node_t &ny) {
    return IndexBlock<node_t>(nx, ny);
  }
  // Construct an instance of near-field block type
  virtual IndexBlock<node_t> makeNearFieldBlock(node_t &nx, node_t &ny) {
    return IndexBlock<node_t>(nx, ny);
  }

 public:
  std::shared_ptr<Tree> rowT;                 // row cluster tree
  std::shared_ptr<Tree> colT;                 // column cluster tree
  std::vector<IndexBlock<node_t>> farField;   // index blocks in the far field
  std::vector<IndexBlock<node_t>> nearField;  // index blocks in the near field
  static bool dbg;                            // Debugging flag
};
/* SAM_LISTING_END_A */

template <class Tree>
bool BlockPartition<Tree>::dbg = false;

/* SAM_LISTING_BEGIN_B */
template <class Tree>
void BlockPartition<Tree>::init(double eta0) {
  buildRec(rowT->root.get(), colT->root.get(), eta0);
}
/* SAM_LISTING_END_B */

/* SAM_LISTING_BEGIN_C */
template <class Tree>
bool BlockPartition<Tree>::adm(const node_t *nx, const node_t *ny,
                               double eta0) const {
  // In an admissible pair neither node must be a leaf.
  if (nx->isLeaf() || ny->isLeaf()) return false;
  // Geometric admissibility condition, see \cref{eq:etadef}.
  const BBox<node_t::dim> Bx = rowT->getBBox(nx), By = colT->getBBox(ny);
  const double bb_dist = dist(Bx, By);
  if (bb_dist == 0.0) return false;
  const double eta = std::max(Bx.diam(), By.diam()) / (2 * bb_dist);
  return (eta < eta0);
}
/* SAM_LISTING_END_C */

/* SAM_LISTING_BEGIN_D */
template <class Tree>
void BlockPartition<Tree>::buildRec(node_t *nx, node_t *ny, double eta0) {
  if (nx && ny) {
    // Add admissible pair to far field
    if (adm(nx, ny, eta0)) {  // \Label[line]{bpr:adm}
      farField.push_back(makeFarFieldBlock(*nx, *ny));
    } else {
      bool rec = false;
      for (int isx = 0; isx <= 1; isx++) {
        for (int isy = 0; isy <= 1; isy++) {
          if (nx->sons[isx] && ny->sons[isy]) {
            // Next level of recursion for non-leaves
            rec = true;
            buildRec(nx->sons[isx].get(), ny->sons[isy].get(), eta0);
          }
        }
      }
      // At least one leaf cluster:
      // Add cluster pair to near field
      if (!rec)  // \Label[line]{bpr:1}
        nearField.push_back(makeNearFieldBlock(*nx, *ny));
    }
  } else
    throw(std::runtime_error("Invalid node pointers"));
}
/* SAM_LISTING_END_D */

// Output near- and far-field block partitioning
template <class Tree>
std::ostream &operator<<(std::ostream &o, const BlockPartition<Tree> &bp) {
  o << "# Near field indices:" << std::endl;
  for (const IndexBlock<typename Tree::node_t> &b : bp.nearField) {
    o << "[ ";
    for (int i : b.i_idx) o << i << " ";
    o << "] x [ ";
    for (int j : b.j_idx) o << j << "  ";
    o << "]" << std::endl;
  }
  o << "# Far field indices:" << std::endl;
  for (const IndexBlock<typename Tree::node_t> &b : bp.farField) {
    o << "[ ";
    for (int i : b.i_idx) o << i << " ";
    o << "] x [";
    for (int j : b.j_idx) o << j << " ";
    o << "]" << std::endl;
  }
  o << "END BLOCK LIST" << std::endl;
  return o;
}

// Auxliary function: Output of matrix block partition in 1D case
void printGeometricPartition(
    const BlockPartition<ClusterTree<CtNode<1>>> &partmat);

}  // namespace HMAT

#endif
