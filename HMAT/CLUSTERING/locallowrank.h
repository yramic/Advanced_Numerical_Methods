/***********************************************************************
 *                                                                     *
 * Code for Course "Advanced Numerical Methods for CSE"                *
 * (Prof. Dr. R. Hiptmair)                                             *
 * Author: Ralf Hiptmair                                               *
 * Date: August 2023                                                   *
 * (C) Seminar for Applied Mathematics, ETH Zurich                     *
 * This code can be freely used for non-commercial purposes as long    *
 * as this header is left intact.                                      *
 ***********************************************************************/

#ifndef LLR_H_
#define LLR_H_

#include <stdexcept>

#include "matrixpartition.h"

namespace HMAT {
// ======================================================================
// For local low-rank approximation
// ======================================================================

/** Node supporting separable approximation by interpolation */
/* SAM_LISTING_BEGIN_Y */
template <int DIM>
class InterpNode : public CtNode<DIM> {
 public:
  // Constructor from sequence of points; initializes \cob{$\VV$}
  InterpNode(const std::vector<Point<DIM>> _pts, std::size_t _q, int _dir = 0)
      : CtNode<DIM>(_pts, _dir), q(_q), sons{nullptr, nullptr} {
    initV();
  }
  virtual ~InterpNode() = default;

 protected:
  // Initialization of matrix \cob{$\VV$}
  void initV();

 public:
  const std::size_t q;               // Rank, no of interpolation nodes
  Eigen::MatrixXd V;                 // low-rank factor \cob{$\VV\in\bbR^{k,q}$}
  std::array<InterpNode *, 2> sons;  // Pointers to sons (of type InterpNode!)
};
/* SAM_LISTING_END_Y */

/* SAM_LISTING_BEGIN_Z */
template <int DIM>
void InterpNode<DIM>::initV() {
  static_assert(DIM == 1, "Implemented only for 1D");
  // **********************************************************************
  // TODO
  // **********************************************************************
}
/* SAM_LISTING_END_Z */

// Recursive output operator for an interpolation node
template <int dim>
std::ostream &operator<<(std::ostream &o, const InterpNode<dim> &node) {
  o << "## IPNode, rank " << node.q << ": " << node.pts.size() << " points, "
    << node.getBBox() << ": ";
  for (const Point<dim> &pt : node.pts) {
    o << "[";
    for (int l = 0; l < dim; l++) {
      o << pt.x[l] << ' ';
    }
    o << "]";
  }
  o << std::endl;
  if (node.sons[0]) o << *(node.sons[0]);
  if (node.sons[1]) o << *(node.sons[1]);
  return o;
}

/** Extended class for cluster trees for local low-rank approximation */
/* SAM_LISTING_BEGIN_E */
template <class Node>
class LLRClusterTree : public ClusterTree<Node> {
 public:
  // Idle constructor just setting rank argument q
  explicit LLRClusterTree(size_t _q) : q(_q) {}
  // Actual constructor taking a sequence of points
  void init(const std::vector<Point<Node::dim>> pts, std::size_t minpts = 1);
  virtual ~LLRClusterTree() = default;

 protected:
  // factory method for relevant type of node taking rank argument
  virtual Node *createNode(const std::vector<Point<Node::dim>> pts, int dir) {
    return new Node(pts, q, dir);
  }

 public:
  const std::size_t q;  // rank of separable approximation on cluster boxes
};

template <class Node>
void LLRClusterTree<Node>::init(const std::vector<Point<Node::dim>> pts,
                                std::size_t minpts) {
  ClusterTree<Node>::init(pts, minpts);
}
/* SAM_LISTING_END_E */

/** Type for far field cluster */
/* SAM_LISTING_BEGIN_F */
template <class Node, typename KERNEL>
class BiDirChebInterpBlock : public IndexBlock<Node> {
 public:
  using kernel_t = KERNEL;
  BiDirChebInterpBlock(const Node &_nx, const Node &_ny, KERNEL _Gfun,
                       std::size_t _q);
  // Constructor that should not be called, needed to avoid compilation errors
  BiDirChebInterpBlock(const Node &_nx, const Node &_ny)
      : IndexBlock<Node>(_nx, _ny), q(0) {
    throw std::runtime_error("Invalid constructor");
  }
  virtual ~BiDirChebInterpBlock() = default;

  KERNEL G;           // kernel function \cob{$\krn$}
  const int q;        // No of interpolation nodes
  Eigen::MatrixXd C;  // \cob{$\VC\in\bbR^{q,q}$}
};
/* SAM_LISTING_END_F */

/* SAM_LISTING_BEGIN_B */
template <class Node, typename KERNEL>
BiDirChebInterpBlock<Node, KERNEL>::BiDirChebInterpBlock(const Node &_nx,
                                                         const Node &_ny,
                                                         KERNEL _Gfun,
                                                         std::size_t _q)
    : IndexBlock<Node>(_nx, _ny), G(std::move(_Gfun)), q(_q) {
  // **********************************************************************
  // TODO
  // **********************************************************************
}
/* SAM_LISTING_END_B */

/** General type for generic near-field cluster pair */
/* SAM_LISTING_BEGIN_Gfun */
template <class Node, typename KERNEL>
class NearFieldBlock : public IndexBlock<Node> {
 public:
  using kernel_t = KERNEL;
  NearFieldBlock(const Node &nx, const Node &ny, KERNEL _Gfun);
  // Constructor that should not be called, needed to avoid compilation errors
  NearFieldBlock(const Node &_nx, const Node &_ny)
      : IndexBlock<Node>(_nx, _ny) {
    throw std::runtime_error("Invalid constructor");
  }

  virtual ~NearFieldBlock() = default;

  KERNEL G;              // kernel function \cob{$\krn$}
  Eigen::MatrixXd Mloc;  // local kernel collocation matrix
};
/* SAM_LISTING_END_Gfun */

/* SAM_LISTING_BEGIN_Q */
template <class Node, typename KERNEL>
NearFieldBlock<Node, KERNEL>::NearFieldBlock(const Node &_nx, const Node &_ny,
                                             KERNEL _Gfun)
    : IndexBlock<Node>(_nx, _ny), G(std::move(_Gfun)) {
  // **********************************************************************
  // TODO
  // **********************************************************************
}
/* SAM_LISTING_END_Q */

/** Extended class for block partition, knowing low-rank compression */
/* SAM_LISTING_BEGIN_H */
template <class Node, typename FFB, typename NFB>
class BiDirChebBlockPartition : public BlockPartition<Node, FFB, NFB> {
 public:
  using kernel_t = typename NFB::kernel_t;
  BiDirChebBlockPartition(std::shared_ptr<const LLRClusterTree<Node>> _xT,
                          std::shared_ptr<const LLRClusterTree<Node>> _yT,
                          kernel_t _Gfun, std::size_t _q, double eta0 = 0.5)
      : BlockPartition<Node, FFB, NFB>(_xT, _yT), G(_Gfun), q(_q) {
    BlockPartition<Node, FFB, NFB>::init(eta0);
  }
  virtual ~BiDirChebBlockPartition() = default;

 protected:
  // Construct an instance of far-field block type
  virtual FFB makeFarFieldBlock(const Node &nx, const Node &ny) const {
    std::cout << "BiDirChebBlockPartition: makeFarFieldBlock" << std::endl;
    return FFB(nx, ny, G, q);
  }
  // Construct an instance of near-field block type
  virtual NFB makeNearFieldBlock(const Node &nx, const Node &ny) const {
    std::cout << "BiDirChebBlockPartition: makeNearFieldBlock" << std::endl;
    return NFB(nx, ny, G);
  }

 public:
  kernel_t G;           // Reference to the kernel function
  const std::size_t q;  // degree+1 of interpolating polynomial
};
/* SAM_LISTING_END_H */

// Special data type for local low-rank compression by one-dimensional
// bi-directional Chebychev interpolation
template <typename KERNEL>
using BiDirChebPartMat1D = HMAT::BiDirChebBlockPartition<
    HMAT::InterpNode<1>, HMAT::BiDirChebInterpBlock<HMAT::CtNode<1>, KERNEL>,
    HMAT::NearFieldBlock<HMAT::CtNode<1>, KERNEL>>;

// Matrix x Vector based on compressed kernel collocation matrix
template <typename KERNEL>
Eigen::VectorXd mvLLRPartMat(const BiDirChebPartMat1D<KERNEL> &Mt,
                             const Eigen::VectorXd &x) {}

}  // namespace HMAT

#endif
