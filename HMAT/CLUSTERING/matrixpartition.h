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

#include <cstddef>
#include <memory>

#include <cassert>

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
  IndexBlock(const NODE &_nx, const NODE &_ny)
      : nx(_nx), ny(_ny), i_idx(_nx.I()), j_idx(ny.I()) {}
  virtual ~IndexBlock() = default;
  const NODE &nx, &ny;                     // contributing clusters
  const std::vector<size_t> i_idx, j_idx;  // contained indices
};
/* SAM_LISTING_END_9 */

/* SAM_LISTING_BEGIN_A */
template <class NODE, typename FFB, typename NFB>
class BlockPartition {
 public:
  using node_t = NODE;
  using farfieldblock_t = FFB;
  using nearfieldblock_t = NFB;
  // Idle constructor
  BlockPartition(std::shared_ptr<const ClusterTree<NODE>> _rowT,
                 std::shared_ptr<const ClusterTree<NODE>> _colT)
      : rowT(_rowT), colT(_colT) {
    assertm((rowT != nullptr), "No valid x-tree!");
    assertm((colT != nullptr), "No valid y-tree!");
  }
  // Trigger recursive construction of partition
  // (Needed, because polymorphic functions not available in constructor)
  void init(double eta0 = 0.5);
  virtual ~BlockPartition() {}
  // Admissibility condition \cob{$\adm$}, see \cref{def:ac}
  virtual bool adm(const NODE *nx, const NODE *ny, double eta0) const;

 protected:
  // Recursive construction from cluster pair
  virtual void buildRec(const NODE *nx, const NODE *ny, double eta0);
  // Construct an instance of far-field block type
  virtual FFB makeFarFieldBlock(const NODE &nx, const NODE &ny) const {
    return FFB(nx, ny);
  }
  // Construct an instance of near-field block type
  virtual NFB makeNearFieldBlock(const NODE &nx, const NODE &ny) const {
    return NFB(nx, ny);
  }

 public:
  std::shared_ptr<const ClusterTree<NODE>> rowT; // row cluster tree
  std::shared_ptr<const ClusterTree<NODE>> colT; // column cluster tree
  std::vector<FFB> farField;   // index blocks in the far field
  std::vector<NFB> nearField;  // index blocks in the near field
  static bool dbg;             // Debugging flag
};
/* SAM_LISTING_END_A */

template <class NODE, typename FFB, typename NFB>
bool BlockPartition<NODE, FFB, NFB>::dbg = false;

/* SAM_LISTING_BEGIN_B */
template <class NODE, typename FFB, typename NFB>
void BlockPartition<NODE, FFB, NFB>::init(double eta0) {
  buildRec(rowT->root, colT->root, eta0);
}
/* SAM_LISTING_END_B */

/* SAM_LISTING_BEGIN_C */
template <class NODE, typename FFB, typename NFB>
bool BlockPartition<NODE, FFB, NFB>::adm(const NODE *nx, const NODE *ny,
                                         double eta0) const {
  // Neither node must be a leaf.
  if (nx->isLeaf() || ny->isLeaf()) return false;
  // Geometric admissibility condition, see \cref{eq:etadef}.
  const BBox<NODE::dim> Bx = nx->getBBox(), By = ny->getBBox();
  const double eta = std::max(Bx.diam(), By.diam()) / (2 * dist(Bx, By));
  return (eta < eta0);
}
/* SAM_LISTING_END_C */

/* SAM_LISTING_BEGIN_D */
template <class NODE, typename FFB, typename NFB>
void BlockPartition<NODE, FFB, NFB>::buildRec(const NODE *nx, const NODE *ny,
                                              double eta0) {
  if (nx && ny) {
    // Add admissible pair to far field
    if (adm(nx, ny, eta0))  // \Label[line]{bpr:adm}
      farField.push_back(makeFarFieldBlock(*nx, *ny));
    else {
      bool rec = false;
      for (int isx = 0; isx <= 1; isx++) {
        for (int isy = 0; isy <= 1; isy++) {
          if (nx->sons[isx] && ny->sons[isy]) {
            // Next level of recursion for non-leaves
            rec = true;
            buildRec(nx->sons[isx], ny->sons[isy], eta0);
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

  // Computation of sparsity measure of a clustertree-based matrix partitioning
  template <typename NODE, typename FFB, typename NFB>
  unsigned int computeSparsityMeasure(const BlockPartition<NODE,FFB,NFB> &bockpart) {
    unsigned int spm = 0;
  // **********************************************************************
  // TODO
  // **********************************************************************
    return spm;
  }
  
template <class Node, typename FFB, typename NFB>
std::ostream &operator<<(std::ostream &o,
                         const BlockPartition<Node, FFB, NFB> &bp) {
  o << "Near field indices:" << std::endl;
  for (const NFB &b : bp.nearField) {
    o << "[ i_idx = ";
    for (int i : b.i_idx) o << i << ", ";
    o << std::endl;
    o << "  j_idx = ";
    for (int j : b.j_idx) o << j << ", ";
    o << "]" << std::endl;
  }
  o << "Far field indices:" << std::endl;
  for (const FFB &b : bp.farField) {
    o << "[ i_idx = ";
    for (int i : b.i_idx) o << i << ", ";
    o << std::endl;
    o << "  j_idx = ";
    for (int j : b.j_idx) o << j << ", ";
    o << "]" << std::endl;
  }
  o << "END BLOCK LIST" << std::endl;
  return o;
}

}  // namespace HMAT

#endif
