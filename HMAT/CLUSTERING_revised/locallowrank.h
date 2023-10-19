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
template <class NODE>
class LLRClusterTree : public ClusterTree<NODE> {
 public:
  // Idle constructor just setting rank argument q
  explicit LLRClusterTree(size_t _q) : q(_q) {}
  // Actual constructor taking a sequence of points
  void init(const std::vector<Point<NODE::dim>> pts, std::size_t minpts = 1);
  virtual ~LLRClusterTree() = default;

 protected:
  // factory method for relevant type of node taking rank argument
  virtual NODE *createNode(const std::vector<Point<NODE::dim>> pts, int dir) {
    return new NODE(pts, q, dir);
  }

 public:
  const std::size_t q;  // rank of separable approximation on cluster boxes
};

template <class NODE>
void LLRClusterTree<NODE>::init(const std::vector<Point<NODE::dim>> pts,
                                std::size_t minpts) {
  ClusterTree<NODE>::init(pts, minpts);
}
/* SAM_LISTING_END_E */

/** Type for far field cluster */
/* SAM_LISTING_BEGIN_F */
template <class NODE, typename KERNEL>
class BiDirChebInterpBlock : public IndexBlock<NODE> {
 public:
  using kernel_t = KERNEL;
  BiDirChebInterpBlock(const NODE &_nx, const NODE &_ny, KERNEL _Gfun,
                       std::size_t _q);
  // Constructor that should not be called, needed to avoid compilation errors
  BiDirChebInterpBlock(const NODE &_nx, const NODE &_ny)
      : IndexBlock<NODE>(_nx, _ny), q(0) {
    throw std::runtime_error("Invalid constructor");
  }
  virtual ~BiDirChebInterpBlock() = default;

  KERNEL G;           // kernel function \cob{$\krn$}
  const int q;        // No of interpolation nodes
  Eigen::MatrixXd C;  // \cob{$\VC\in\bbR^{q,q}$}
};
/* SAM_LISTING_END_F */

/* SAM_LISTING_BEGIN_B */
template <class NODE, typename KERNEL>
BiDirChebInterpBlock<NODE, KERNEL>::BiDirChebInterpBlock(const NODE &_nx,
                                                         const NODE &_ny,
                                                         KERNEL _Gfun,
                                                         std::size_t _q)
    : IndexBlock<NODE>(_nx, _ny), G(std::move(_Gfun)), q(_q) {
  // **********************************************************************
  // TODO
  // **********************************************************************
}
/* SAM_LISTING_END_B */

/** General type for generic near-field cluster pair */
/* SAM_LISTING_BEGIN_G */
template <class NODE, typename KERNEL>
class NearFieldBlock : public IndexBlock<NODE> {
 public:
  using kernel_t = KERNEL;
  NearFieldBlock(const NODE &nx, const NODE &ny, KERNEL _Gfun);
  // Constructor that should not be called, needed to avoid compilation errors
  NearFieldBlock(const NODE &_nx, const NODE &_ny)
      : IndexBlock<NODE>(_nx, _ny) {
    throw std::runtime_error("Invalid constructor");
  }

  virtual ~NearFieldBlock() = default;

  KERNEL G;              // kernel function \cob{$\krn$}
  Eigen::MatrixXd Mloc;  // local kernel collocation matrix
};
/* SAM_LISTING_END_G */

/* SAM_LISTING_BEGIN_Q */
template <class NODE, typename KERNEL>
NearFieldBlock<NODE, KERNEL>::NearFieldBlock(const NODE &_nx, const NODE &_ny,
                                             KERNEL _Gfun)
    : IndexBlock<NODE>(_nx, _ny), G(std::move(_Gfun)) {
  // **********************************************************************
  // TODO
  // **********************************************************************
}
/* SAM_LISTING_END_Q */

/** Extended class for block partition, knowing low-rank compression */
/* SAM_LISTING_BEGIN_H */
template <class NODE, typename FFB, typename NFB>
class BiDirChebBlockPartition : public BlockPartition<NODE, FFB, NFB> {
 public:
  using kernel_t = typename NFB::kernel_t;
  BiDirChebBlockPartition(std::shared_ptr<const LLRClusterTree<NODE>> _xT,
                          std::shared_ptr<const LLRClusterTree<NODE>> _yT,
                          kernel_t _Gfun, std::size_t _q, double eta0 = 0.5)
      : BlockPartition<NODE, FFB, NFB>(_xT, _yT), G(_Gfun), q(_q) {
    BlockPartition<NODE, FFB, NFB>::init(eta0);
  }
  virtual ~BiDirChebBlockPartition() = default;

 protected:
  // Construct an instance of far-field block type
  virtual FFB makeFarFieldBlock(const NODE &nx, const NODE &ny) const {
    std::cout << "BiDirChebBlockPartition: makeFarFieldBlock" << std::endl;
    return FFB(nx, ny, G, q);
  }
  // Construct an instance of near-field block type
  virtual NFB makeNearFieldBlock(const NODE &nx, const NODE &ny) const {
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
/* SAM_LISTING_BEGIN_U */
template <typename PARTMAT1D>
Eigen::VectorXd mvLLRPartMat(const PARTMAT1D &Mt, const Eigen::VectorXd &x) {
  // **********************************************************************
  // TODO
  // **********************************************************************
}
/* SAM_LISTING_END_U */

// Computing scaled Frobenius norm for local low-rank matrix approximation error
/* SAM_LISTING_BEGIN_V */
template <typename KERNEL>
double approxErrorLLR(const BiDirChebPartMat1D<KERNEL> &Mt) {
  // **********************************************************************
  // TODO
  // *********************************************************************
}
/* SAM_LISTING_END_V */

// Validation of implementation of local low-rank compression based on
// bi-directional Chebychev interpolation
bool validateLLR(unsigned int q);

// Tabulate approximation errors
void tabulateConvergenceLLR(std::vector<unsigned int> &&n,
                            std::vector<unsigned int> &&q);

// Measure runtimes
void runtimeMatVec(std::vector<unsigned int> &&n);

/* SAM_LISTING_BEGIN_R */
template <int d>  // dimension \cob{$d$} as template argument
struct BasisFn {
  std::size_t idx;                         // Index of basis function
  Eigen::Matrix<double, d, 1> xmin, xmax;  // Corners of bounding box
};
/* SAM_LISTING_END_R */

}  // namespace HMAT

#endif
