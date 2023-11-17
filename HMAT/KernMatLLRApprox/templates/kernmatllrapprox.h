/**
 * @file kernmatllrapprox.h
 * @brief NPDE homework KernMatLLRApprox code
 * @author R. Hiptmair
 * @date September 2023
 * @copyright Developed at SAM, ETH Zurich
 */

#ifndef KERNMATLLRAPPROX_H_
#define KERNMATLLRAPPROX_H_
#define _USE_MATH_DEFINES
#include <matrixpartition.h>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <iomanip>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <unordered_map>

namespace KernMatLLRApprox {

/** @brief verify that a ClusterTree object complies with the definition of a
cluster tree
*
*/
/* SAM_LISTING_BEGIN_1 */
// clang-format off
template <typename NODE>
bool checkClusterTree(const HMAT::ClusterTree<NODE> &T) {
  if (T.root == nullptr) return false;
  bool ok = true;  // Flag to be returned
  // **********************************************************************
  // TO BE SUPPLEMENTED
  // **********************************************************************
  return ok;
}
/* SAM_LISTING_END_1 */
// clang-format on

// clang-format off
/** @brief Verify validity of matrix paving into near and far field blocks */
/* SAM_LISTING_BEGIN_2 */
template <typename NODE>
bool checkMatrixPartition(
    const HMAT::BlockPartition<NODE, HMAT::IndexBlock<NODE>,
                               HMAT::IndexBlock<NODE>> &partmat) {
  assertm((partmat.rowT and partmat.colT), "Missing trees!");
  // **********************************************************************
  // TO BE SUPPLEMENTED
  // **********************************************************************
  return false;
}
/* SAM_LISTING_END_2 */
// clang-format on

// ======================================================================
// For local low-rank approximation
// ======================================================================

/** Node supporting separable approximation by interpolation and also efficient
   matrix x vector operations */
// clang-format off
/* SAM_LISTING_BEGIN_Y */
template <int DIM>
class InterpNode : public HMAT::CtNode<DIM> {
 public:
  // Constructor from sequence of points; initializes \cob{$\VV$}
  InterpNode(const std::vector<HMAT::Point<DIM>> _pts, std::size_t _q,
             int _dir = 0)
      : HMAT::CtNode<DIM>(_pts, _dir), q(_q),  sons{nullptr, nullptr},
    clust_omega(_q), V(_pts.size(), _q) { initV(); }
  virtual ~InterpNode() = default;
  // Is the node a leaf node ?
  [[nodiscard]] bool isLeaf() const override {
    return (!(sons[0]) and !(sons[1]));
  }
 protected:
  // Initialization of the low-rank factor matrix \cob{$\VV_w$}
  void initV();
 public:
  const std::size_t q;               // Rank, no of interpolation nodes
  Eigen::MatrixXd V;                 // low-rank factor \cob{$\VV\in\bbR^{k,q}$}
  std::array<InterpNode *, 2> sons;  // Pointers to sons (of type InterpNode!)
  Eigen::VectorXd clust_omega;     // for cluster-local linear algebra
};
/* SAM_LISTING_END_Y */
// clang-format on

// Recursive output operator for an interpolation node
template <int dim>
std::ostream &operator<<(std::ostream &o, const InterpNode<dim> &node) {
  o << "## IPNode, rank " << node.q << ": " << node.pts.size() << " points, ";
  if (node.sons[0]) o << "S0 ";
  if (node.sons[1]) o << "S1 ";
  if (!(node.sons[0]) && !(node.sons[1])) {
    o << "LEAF  ";
  }
  o << node.getBBox() << ": ";
  for (const HMAT::Point<dim> &pt : node.pts) {
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

// For debugging purposes: non-recursive output for 1D
void outIPNode(const InterpNode<1> &ipnode, bool printV = false,
               std::ostream &o = std::cout);
// clang-format off
/* SAM_LISTING_BEGIN_Z */
template <int DIM>
void InterpNode<DIM>::initV() {
  static_assert(DIM == 1, "Implemented only for 1D");
  // Retrieve bounding box, \href{https://stackoverflow.com/questions/4643074/why-do-i-have-to-access-template-base-class-members-through-the-this-pointer}{eplanation} for the this pointer
  const HMAT::BBox<DIM> bbox(this->pts);
  // Find interval correspoding to the bounding box of the current cluster
  const double a = bbox.minc[0];
  const double b = bbox.maxc[0];
    // **********************************************************************
    // TODO
    // **********************************************************************
}
/* SAM_LISTING_END_Z */
// clang-format on

/** Extended class for cluster trees for local low-rank approximation */
/* SAM_LISTING_BEGIN_E */
template <class NODE>
class LLRClusterTree : public HMAT::ClusterTree<NODE> {
 public:
  // Idle constructor just setting rank argument q
  explicit LLRClusterTree(size_t _q) : q(_q) {}
  // Actual constructor taking a sequence of points
  void init(const std::vector<HMAT::Point<NODE::dim>> pts,
            std::size_t minpts = 1);
  virtual ~LLRClusterTree() = default;

 protected:
  // factory method for relevant type of node taking rank argument
  virtual NODE *createNode(const std::vector<HMAT::Point<NODE::dim>> pts,
                           int dir) {
    return new NODE(pts, q, dir);
  }

 public:
  const std::size_t q;  // rank of separable approximation on cluster boxes
};

template <class NODE>
void LLRClusterTree<NODE>::init(const std::vector<HMAT::Point<NODE::dim>> pts,
                                std::size_t minpts) {
  HMAT::ClusterTree<NODE>::init(pts, minpts);
}
/* SAM_LISTING_END_E */

/** Type for far field cluster */
/* SAM_LISTING_BEGIN_F */
template <class NODE, typename KERNEL>
class BiDirChebInterpBlock : public HMAT::IndexBlock<NODE> {
 public:
  using kernel_t = KERNEL;
  BiDirChebInterpBlock(NODE &_nx, NODE &_ny, KERNEL _Gfun, std::size_t _q);
  // Constructor that should not be called, needed to avoid compilation errors
  BiDirChebInterpBlock(NODE &_nx, NODE &_ny)
      : HMAT::IndexBlock<NODE>(_nx, _ny), q(0) {
    throw std::runtime_error("Invalid constructor");
  }
  virtual ~BiDirChebInterpBlock() = default;

  KERNEL G;           // kernel function \cob{$\krn$}
  const int q;        // No of interpolation nodes
  Eigen::MatrixXd C;  // \cob{$\VC\in\bbR^{q,q}$}, see \lref{eq:rqbd}
};
/* SAM_LISTING_END_F */

// clang-format off
/* SAM_LISTING_BEGIN_B */
template <class NODE, typename KERNEL>
BiDirChebInterpBlock<NODE, KERNEL>::BiDirChebInterpBlock(
 NODE &_nx, NODE &_ny, KERNEL _Gfun, std::size_t _q)
    : HMAT::IndexBlock<NODE>(_nx, _ny), G(std::move(_Gfun)), q(_q), C(_q, _q) {
  static_assert(NODE::dim == 1, "Only implemented in 1D");
    // **********************************************************************
    // TODO
    // **********************************************************************
}
/* SAM_LISTING_END_B */
// clang-format on

/** General type for generic near-field cluster pair */
/* SAM_LISTING_BEGIN_G */
template <class NODE, typename KERNEL>
class NearFieldBlock : public HMAT::IndexBlock<NODE> {
 public:
  using kernel_t = KERNEL;
  NearFieldBlock(NODE &nx, NODE &ny, KERNEL _Gfun);
  // Constructor that should not be called, needed to avoid compilation errors
  NearFieldBlock(NODE &_nx, NODE &_ny) : HMAT::IndexBlock<NODE>(_nx, _ny) {
    throw std::runtime_error("Invalid constructor");
  }

  virtual ~NearFieldBlock() = default;

  KERNEL G;              // kernel function \cob{$\krn$}
  Eigen::MatrixXd Mloc;  // local kernel collocation matrix
};
/* SAM_LISTING_END_G */

/* SAM_LISTING_BEGIN_Q */
template <class NODE, typename KERNEL>
NearFieldBlock<NODE, KERNEL>::NearFieldBlock(NODE &_nx, NODE &_ny, KERNEL _Gfun)
    : HMAT::IndexBlock<NODE>(_nx, _ny),
      G(std::move(_Gfun)),
      Mloc(_nx.pts.size(), _ny.pts.size()) {
  static_assert(NODE::dim == 1, "Only implemented in 1D");
  // **********************************************************************
  // TODO
  // **********************************************************************
}
/* SAM_LISTING_END_Q */

/** Extended class for block partition, knowing low-rank compression */
/* SAM_LISTING_BEGIN_H */
template <class NODE, typename FFB, typename NFB>
class BiDirChebBlockPartition : public HMAT::BlockPartition<NODE, FFB, NFB> {
 public:
  using kernel_t = typename NFB::kernel_t;
  BiDirChebBlockPartition(std::shared_ptr<LLRClusterTree<NODE>> _rowT,
                          std::shared_ptr<LLRClusterTree<NODE>> _colT,
                          kernel_t _Gfun, std::size_t _q, double eta0 = 2.0)
      : HMAT::BlockPartition<NODE, FFB, NFB>(_rowT, _colT), G(_Gfun), q(_q) {
    HMAT::BlockPartition<NODE, FFB, NFB>::init(eta0);
  }
  virtual ~BiDirChebBlockPartition() = default;

 protected:
  // Construct an instance of far-field block type
  virtual FFB makeFarFieldBlock(NODE &nx, NODE &ny) {
    ffb_cnt++;
    return FFB(nx, ny, G, q);
  }
  // Construct an instance of near-field block type
  virtual NFB makeNearFieldBlock(NODE &nx, NODE &ny) {
    nfb_cnt++;
    return NFB(nx, ny, G);
  }

 public:
  kernel_t G;               // Reference to the kernel function
  const std::size_t q;      // degree+1 of interpolating polynomial
  unsigned int ffb_cnt{0};  // Counter for far-field blocks
  unsigned int nfb_cnt{0};  // Counter for near-field blocks
};
/* SAM_LISTING_END_H */

// Special data type for local low-rank compression by one-dimensional
// bi-directional Chebychev interpolation
template <typename KERNEL>
using BiDirChebPartMat1D =
    BiDirChebBlockPartition<InterpNode<1>,
                            BiDirChebInterpBlock<InterpNode<1>, KERNEL>,
                            NearFieldBlock<HMAT::CtNode<1>, KERNEL>>;

// clang-format off
// Matrix x Vector based on compressed kernel collocation matrix
// Assumes contiguous indices in root node!
// The NODE template parameter must be compatible with InterpNode
/* SAM_LISTING_BEGIN_U */
template <class NODE, typename FFB, typename NFB>
Eigen::VectorXd mvLLRPartMat(BiDirChebBlockPartition<NODE, FFB, NFB> &llrcmat,
                             const Eigen::VectorXd &x) {
  using ff_node_t = typename FFB::node_t;
  using nf_node_t = typename NFB::node_t;
  // Verify requirement of contiguous indices
  const std::vector<size_t> col_idxs = (llrcmat.colT->root)->I();
  assertm((*std::max_element(col_idxs.begin(), col_idxs.end()) ==
           col_idxs.size() - 1),
          "col idxs not contiguous");
  const std::vector<size_t> row_idxs = (llrcmat.rowT->root)->I();
  assertm((*std::max_element(row_idxs.begin(), row_idxs.end()) ==
           row_idxs.size() - 1),
          "row idxs not contiguous");
  assertm(x.size() == col_idxs.size(), "Wrong size of vector x");
  Eigen::VectorXd y(row_idxs.size());  // Return value
  y.setZero();
// **********************************************************************
// TODO
// **********************************************************************
  return y;
}
/* SAM_LISTING_END_U */
// clang-format on

// Computing scaled Frobenius norm for local low-rank matrix approximation error
/* SAM_LISTING_BEGIN_V */
template <typename KERNEL>
std::pair<double, double> approxErrorLLR(BiDirChebPartMat1D<KERNEL> &llrcM) {
  const size_t n = llrcM.rows();
  const size_t m = llrcM.cols();
  Eigen::MatrixXd M(n, m);   // Exact kernel collocation matrix
  Eigen::MatrixXd Mt(n, m);  // Compressed matrix as dense matrix
      // **********************************************************************
      // TODO
      // *********************************************************************
  return {std::sqrt((M - Mt).squaredNorm() / (n * m)),
          std::sqrt(M.squaredNorm() / (n * m))};
}
/* SAM_LISTING_END_V */

/* SAM_LISTING_BEGIN_S */
template <typename NODE, typename FFB, typename NFB>
unsigned int computeSparsityMeasure(
    const HMAT::BlockPartition<NODE, FFB, NFB> &blockpart,
    std::ostream *out = nullptr) {
  assertm((blockpart.rowT and blockpart.colT), "Missing trees!");
  // Set up hash maps for nodes of both trees
  using nf_node_t = typename NFB::node_t;
  using hashmap_t = std::unordered_map<const nf_node_t *, int>;
  using keyval_t = typename hashmap_t::value_type;
  hashmap_t nodemap_row;
  hashmap_t nodemap_col;
  // Maximal node counts for row adn column clusters
  int xnode_maxcnt = 0;
  int ynode_maxcnt = 0;
  // **********************************************************************
  // YOUR CODE HERE
  // **********************************************************************
  return std::max(xnode_maxcnt, ynode_maxcnt);
}
/* SAM_LISTING_END_S */

// Validation of implementation of local low-rank compression based on
// bi-directional Chebychev interpolation
bool validateLLR(unsigned int q, double tol = 1.0E-8, double eta = 2.0);

// Tabulate approximation errors
void tabulateConvergenceLLR(std::vector<unsigned int> &&n_vec,
                            std::vector<unsigned int> &&q_vec,
                            double eta = 2.0);

// Measure runtimes
void runtimeMatVec(std::vector<unsigned int> &&n_vec, unsigned int n_runs = 3,
                   unsigned int q = 5, double eta = 2.0);

}  // namespace KernMatLLRApprox

#endif
