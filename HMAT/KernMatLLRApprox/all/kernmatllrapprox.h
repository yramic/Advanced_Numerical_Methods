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
#include <stdexcept>

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
#if SOLUTION
  // Index set $\mathbb{I}$ underlying the cluster tree
  const std::vector<std::size_t> idxset{T.root->I()};
  // Array for remapping indices to integers from 0 to $\sharp\mathbb{I}$
  const std::size_t maxidx = *std::max_element(idxset.begin(), idxset.end());
  // \cor{A bit of a gamble: indices must not be too large compared to their number}
  assertm((maxidx < 2 * idxset.size()), "Maximal index too large");
  std::vector<int> remap_idx(maxidx + 1, -1);
  for (int k = 0; k < idxset.size(); ++k) remap_idx[idxset[k]] = k;
  std::vector<unsigned int> c(idxset.size(), 0);  // Work array

  // Recursive checking
  std::function<void(const NODE *)> ctcheck_rec =
      [&](const NODE *node) -> void {
    if (node == nullptr) return;
    // Indices present in the parent node are marked in the work array
    std::fill(c.begin(), c.end(), 0);  // Initialize work array with zero
    const std::vector<std::size_t> idxset_node{node->I()};  // Index set of node
    for (std::size_t idx : idxset_node) {
      if (idx > maxidx) { ok = false;  return; } // Index not present in root
      else {
        const int pos = remap_idx[idx];
        if (pos < 0) { ok = false; return; } // index not in root index set
        c[pos] = 1; // mark index as present in parent node
      }
    }
    // Visit all sons
    bool has_sons = false; // no sons at all ? 
    for (const NODE *son : node->sons) {
      if (son != nullptr) {
        has_sons = true;
        const std::vector<std::size_t> idxset_son{son->I()}; // son's index set 
        // Increment entries of work array for indices present in a son
        for (std::size_t idx : idxset_son) {
          if (idx > maxidx) { ok = false; return; } // Index not present in root 
	  else {
            const int pos = remap_idx[idx];
            if (pos < 0) { ok = false; return; } // Index not present in root 
            c[pos] += 1; 
          } } } }
    // The work array may only contain zeros or 2s at this point
    if (has_sons) {
      for (unsigned int c_ent : c) {
        if ((c_ent != 0) and (c_ent != 2)) { ok = false; return; }
      }
      // Recursive call
      if (ok and has_sons) {
        for (const NODE *son : node->sons) {
	  ctcheck_rec(son);
	  if (!ok) return;
	} } }
  };
  // Lanch recursion
  ctcheck_rec(T.root);
#else
  // **********************************************************************
  // TO BE SUPPLEMENTED
  // **********************************************************************
#endif
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
#if SOLUTION
  // Index sets $\mathbb{I}$, $\mathbb{J}$ underlying the row/col cluster trees
  const std::vector<std::size_t> idxset_row{partmat.rowT->root->I()};
  const std::vector<std::size_t> idxset_col{partmat.colT->root->I()};
  // Array for remapping indices to integers from 0 to
  // $\sharp\mathbb{I}/\mathbb{J}$
  const std::size_t maxidx_row =
      *std::max_element(idxset_row.begin(), idxset_row.end());
  const std::size_t maxidx_col =
      *std::max_element(idxset_col.begin(), idxset_col.end());
  // \cor{A bit of a gamble: indices $\approx\sharp\mathbb{I},\sharp\mathbb{J}$}
  assertm((maxidx_row < 2 * idxset_row.size()) and
              (maxidx_col < 2 * idxset_col.size()),
          "Maximal index too large");
  std::vector<int> remap_idx_row(maxidx_row + 1, -1);
  std::vector<int> remap_idx_col(maxidx_col + 1, -1);
  for (int k = 0; k < idxset_row.size(); ++k) remap_idx_row[idxset_row[k]] = k;
  for (int k = 0; k < idxset_col.size(); ++k) remap_idx_col[idxset_col[k]] = k;
  // Matrix for keeping track of occurrences of matrix entries
  Eigen::MatrixXi C(idxset_row.size(), idxset_col.size());
  C.setZero();

  // Run through all far-field blocks
  for (const auto &ffb : partmat.farField) {
    for (std::size_t row_idx : ffb.i_idx) {
      for (std::size_t col_idx : ffb.j_idx) {
        if ((row_idx > maxidx_row) or (col_idx > maxidx_col)) return false;
        const int pos_row = remap_idx_row[row_idx];
        const int pos_col = remap_idx_col[col_idx];
        if ((pos_row < 0) or (pos_col < 0)) return false;
        C(pos_row, pos_col) += 1;
      }
    }
  }
  // Run through all near-field blocks: samne code
  for (const auto &nfb : partmat.nearField) {
    for (std::size_t row_idx : nfb.i_idx) {
      for (std::size_t col_idx : nfb.j_idx) {
        if ((row_idx > maxidx_row) or (col_idx > maxidx_col)) return false;
        const int pos_row = remap_idx_row[row_idx];
        const int pos_col = remap_idx_col[col_idx];
        if ((pos_row < 0) or (pos_col < 0)) return false;
        C(pos_row, pos_col) += 1;
      }
    }
  }
  // Check whether matrix contains entries != 1
  return (C.array() == 1).all();
#else
  // **********************************************************************
  // TO BE SUPPLEMENTED
  // **********************************************************************
  return false;
#endif
}
/* SAM_LISTING_END_2 */
// clang-format on

// ======================================================================
// For local low-rank approximation
// ======================================================================

/** Node supporting separable approximation by interpolation */
/* SAM_LISTING_BEGIN_Y */
template <int DIM>
class InterpNode : public HMAT::CtNode<DIM> {
 public:
  // Constructor from sequence of points; initializes \cob{$\VV$}
  InterpNode(const std::vector<HMAT::Point<DIM>> _pts, std::size_t _q,
             int _dir = 0)
      : HMAT::CtNode<DIM>(_pts, _dir), q(_q), sons{nullptr, nullptr} {
    initV();
  }
  virtual ~InterpNode() = default;
  // protected:
  // Initialization of the low-rank factor matrix \cob{$\VV_w$}
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
  // Retrieve bounding box,
  // \href{https://stackoverflow.com/questions/4643074/why-do-i-have-to-access-template-base-class-members-through-the-this-pointer}{eplanation}
  // for the this pointer
  const HMAT::BBox<DIM> bbox(this->pts);
  // Find interval correspoding to the bounding box of the current cluster
  const double a = bbox.minc[0];
  const double b = bbox.maxc[0];
#if SOLUTION
  // Compute Chebychev nodes $t_i$ and
  // baryccentric weights $\lambda_i$ for interval $\cintv{a,b}$
  Eigen::VectorXd lambda(q);  // barycentric weights
  Eigen::VectorXd t(q);       // Chebychev nodes
  const double fac = std::pow(4.0 / (b - a), q - 1) / q;
  int sgn = 1;
  for (int i = 0; i < q; i++, sgn *= -1) {
    const double arg = (2.0 * i + 1.0) / (2 * q) * M_PI;
    t[i] = a + 0.5 * (b - a) * (std::cos(arg) + 1.0);
    lambda[i] = fac * sgn * std::sin(arg);
  }
  // Number of collocation points in cluster
  const unsigned int nIw = (this->pts).size();
  Eigen::MatrixXd V(nIw, q);
  // Traverse collocation points (in the interval $\cintv{a,b}$)
  Eigen::VectorXd tx_diff(q);  // $t_i - \xi$
  double tau = 0.0;  // $\cob{tau := sum_{i=1}^{q} \frac{\lambda_i}{\xi-t_i}}$
  for (int j = 0; j < nIw; ++j) {
    // Coordinate of current collocation point
    const double xi = ((this->pts)[j]).x[0];
    bool on_node = false;
    for (int i = 0; i < q; i++) {
      tx_diff[i] = xi - t[i];
      // Avoid division by zero
      if (tx_diff[i] == 0.0) {
        on_node = true;
        V.row(j).setZero();
        V(j, i) = 1.0;  // $\cob{(\VV)_{j,:} = \Ve_{i}^{\top}}$
        break;          // The $j$-th row of $\VV$ is complete already
      } else {
        const double txdl = lambda[i] / tx_diff[i];
        V(j, i) = txdl;
        tau += txdl;
      }
    }
    if (!on_node) {
      V.row(j) /= tau;
    }
  }
#else
// **********************************************************************
// TODO
// **********************************************************************
#endif
}
/* SAM_LISTING_END_Z */

// Recursive output operator for an interpolation node
template <int dim>
std::ostream &operator<<(std::ostream &o, const InterpNode<dim> &node) {
  o << "## IPNode, rank " << node.q << ": " << node.pts.size() << " points, "
    << node.getBBox() << ": ";
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
  BiDirChebInterpBlock(const NODE &_nx, const NODE &_ny, KERNEL _Gfun,
                       std::size_t _q);
  // Constructor that should not be called, needed to avoid compilation errors
  BiDirChebInterpBlock(const NODE &_nx, const NODE &_ny)
      : HMAT::IndexBlock<NODE>(_nx, _ny), q(0) {
    throw std::runtime_error("Invalid constructor");
  }
  virtual ~BiDirChebInterpBlock() = default;

  KERNEL G;           // kernel function \cob{$\krn$}
  const int q;        // No of interpolation nodes
  Eigen::MatrixXd C;  // \cob{$\VC\in\bbR^{q,q}$}, see \lref{eq:rqbd}
};
/* SAM_LISTING_END_F */

/* SAM_LISTING_BEGIN_B */
template <class NODE, typename KERNEL>
BiDirChebInterpBlock<NODE, KERNEL>::BiDirChebInterpBlock(const NODE &_nx,
                                                         const NODE &_ny,
                                                         KERNEL _Gfun,
                                                         std::size_t _q)
    : HMAT::IndexBlock<NODE>(_nx, _ny), G(std::move(_Gfun)), q(_q) {
  static_assert(NODE::dim == 1, "Only implemented in 1D");
  // Obtain bounding boxes, here intervals
  std::array<HMAT::BBox<1>, 2> bboxs = {_nx.getBBox(), _nx.getBBox()};
  std::array<Eigen::Vector2d, 2> intv;
  for (int d = 0; d < 2; ++d) {
    intv[d] = Eigen::Vector2d(bboxs[d].minc[0], bboxs[d].maxc[0]);
  }
  // Compute Chebychev nodes
  Eigen::MatrixXd t(2, q);
  for (int j = 0; j < q; ++j) {
    const double cosval = std::cos((2.0 * j + 1.0) / (2 * q) * M_PI);
    for (int d = 0; d < 2; ++d) {
      t(d, j) = intv[d][0] + 0.5 * (intv[d][1] - intv[d][0]) * (cosval + 1.0);
    }
  }
  // Fill matrix $\VC$ 

  // **********************************************************************
  // TODO
  // **********************************************************************
}
/* SAM_LISTING_END_B */

/** General type for generic near-field cluster pair */
/* SAM_LISTING_BEGIN_G */
template <class NODE, typename KERNEL>
class NearFieldBlock : public HMAT::IndexBlock<NODE> {
 public:
  using kernel_t = KERNEL;
  NearFieldBlock(const NODE &nx, const NODE &ny, KERNEL _Gfun);
  // Constructor that should not be called, needed to avoid compilation errors
  NearFieldBlock(const NODE &_nx, const NODE &_ny)
      : HMAT::IndexBlock<NODE>(_nx, _ny) {
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
    : HMAT::IndexBlock<NODE>(_nx, _ny), G(std::move(_Gfun)) {
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
  BiDirChebBlockPartition(std::shared_ptr<const LLRClusterTree<NODE>> _rowT,
                          std::shared_ptr<const LLRClusterTree<NODE>> _colT,
                          kernel_t _Gfun, std::size_t _q, double eta0 = 0.5)
      : HMAT::BlockPartition<NODE, FFB, NFB>(_rowT, _colT), G(_Gfun), q(_q) {
    HMAT::BlockPartition<NODE, FFB, NFB>::init(eta0);
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
using BiDirChebPartMat1D =
    BiDirChebBlockPartition<InterpNode<1>,
                            BiDirChebInterpBlock<HMAT::CtNode<1>, KERNEL>,
                            NearFieldBlock<HMAT::CtNode<1>, KERNEL>>;

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

}  // namespace KernMatLLRApprox

#endif
