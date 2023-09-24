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
  // Compute Chebychev nodes $t_i$ and
  // baryccentric weights $\lambda_i$ for interval $\cintv{a,b}$
  Eigen::VectorXd lambda(q);  // barycentric weights
  Eigen::VectorXd t(q);       // Chebychev nodes
  const double fac = std::pow(4.0 / (b - a), q - 1) / q;
  int sgn = 1;
  for (int i = 0; i < q; i++, sgn *= -1) {
    const double arg = (2.0 * i + 1.0) / (2 * q) * M_PI;
    t[i] = a + 0.5 * (b - a) * (std::cos(arg) + 1.0);
    lambda[i] = fac * sgn * std::sin(arg); // \prbeqref{eq:bwf}
  }
  // Number of collocation points in cluster
  const unsigned int nIw = (this->pts).size();
  // Traverse collocation points (in the interval $\cintv{a,b}$)
  Eigen::VectorXd tx_diff(q);  // $t_i - \xi$
  for (int j = 0; j < nIw; ++j) {
    // Coordinate of current collocation point
    const double xi = ((this->pts)[j]).x[0];
    double tau = 0.0;  // $\cob{\tau := \sum_{i=1}^{q} \frac{\lambda_i}{\xi-t_i}}$
    bool on_node = false;
    for (int i = 0; i < q; i++) {
      tx_diff[i] = xi - t[i];
      // Avoid division by zero
      if (tx_diff[i] == 0.0) {
        on_node = true;
        V.row(j).setZero();
        V(j, i) = 1.0;  // $\cob{(\VV)_{j,:} = \Ve_{i}^{\top}}$ when hitting a node
        break;          // The $j$-th row of $\VV$ is complete already
      } else {
        const double txdl = lambda[i] / tx_diff[i];
        V(j, i) = txdl;
        tau += txdl;
      }
    }
    if (!on_node) V.row(j) /= tau;
  }
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
  // Obtain bounding boxes, here intervals
  std::array<HMAT::BBox<1>, 2> bboxs = {_nx.getBBox(), _ny.getBBox()};
  std::array<Eigen::Vector2d, 2> intv;
  for (int d = 0; d < 2; ++d) {
    intv[d] = Eigen::Vector2d(bboxs[d].minc[0], bboxs[d].maxc[0]);
  }
  // Compute Chebychev nodes for nodal intervals, \lref{eq:chn}
  Eigen::MatrixXd t(2, q);
  for (int j = 0; j < q; ++j) {
    // Chebychev nodes on standard interval 
    const double cosval = std::cos((2.0 * j + 1.0) / (2 * q) * M_PI);
    for (int d = 0; d < 2; ++d) {
      // Code for \prbeqref{eq:tvk} \& \prbeqref{eq:twl}
      t(d, j) = intv[d][0] + 0.5 * (intv[d][1] - intv[d][0]) * (cosval + 1.0);
    }
  }
  // Fill matrix $\cob{\VC}$: $\cob{(\VC)_{k,\ell} = G(t_k,t_{\ell})}$, see \prbeqref{eq:Cm}
  for (int k = 0; k < q; ++k) {
    for (int j = 0; j < q; ++j) {
      C(k, j) = G(t(0, k), t(1, j));
    }
  }
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
// Direct initialization of near field kernel collocation matrix
  for (int i = 0; i < _nx.pts.size(); ++i) {
    for (int j = 0; j < _ny.pts.size(); ++j) {
      Mloc(i, j) = G((_nx.pts[i]).x[0], (_ny.pts[j]).x[0]);
    }
  }
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
  // Pass (I) of \lref{par:3p}: Reduce to cluster for column tree, computation of
  // $\cob{\vec{\omegabf}_{w} := \cog{\VV_{w}^{\top}}\VR_{w}\vec{\mubf}\in\bbR^{\sharp{\Ci(w)}}}$
  std::function<void(NODE * col_node)> comp_omega_rec =
      [&](NODE *col_node) -> void {
    if (col_node) {
      // Obtain indices held by the cluster
      const std::vector<size_t> idxs = col_node->I();
      for (int j = 0; j < idxs.size(); ++j) {
        // Restriction of argument vector to column cluster
        (col_node->clust_sect_vec)[j] = x[idxs[j]];
        // Compute $\cob{\vec{\omegabf}_{w}}$
        col_node->clust_omega =
            (col_node->V).transpose() * col_node->clust_sect_vec;
      }
      // Recursion: visit entire cluster tree
      comp_omega_rec(col_node->sons[0]);
      comp_omega_rec(col_node->sons[1]);
    }
  };
  comp_omega_rec(llrcmat.colT->root);

  // Final sentence in \lref{par:3p}: Clear local storage of row tree
  std::function<void(NODE * ipnode)> clear_vec_rec = [&](NODE *ipnode) -> void {
    if (ipnode) {
      ipnode->clust_sect_vec.setZero();
      ipnode->clust_omega.setZero();
      clear_vec_rec(ipnode->sons[0]);
      clear_vec_rec(ipnode->sons[1]);
    }
  };
  clear_vec_rec(llrcmat.rowT->root);

  // Pass (II) in \lref{par:3p}: Do block-based computations
  // Visit far-field blocks:
  for (auto &ffb : llrcmat.farField) {
    // Far-field block: accumulate $\cob{\VC_{v\times w}\vec{\omegabf}_w}$
    ff_node_t &row_node = ffb.nx;
    const ff_node_t &col_node = ffb.ny;
    row_node.clust_omega += ffb.C * col_node.clust_omega;
  }
  for (auto &nfb : llrcmat.nearField) {
    // Near-field block: use retricted kernel collocatin matrix $\cob{\rst{\VM}{v\times w}}$
    nf_node_t &row_node = nfb.nx;
    const nf_node_t &col_node = nfb.ny;
    row_node.clust_sect_vec += nfb.Mloc * col_node.clust_sect_vec;
  }

  // Pass (III) of \lref{par:3p}: Traverse row tree and assemble return vector
  std::function<void(NODE * row_node)> ass_y_rec = [&](NODE *row_node) -> void {
    if (row_node) {
      const std::vector<size_t> idxs = row_node->I();
      // $\cob{\rst{\Vx}{v} += \VU_{v}\vec{\zetabf}_{v}+\vec{\phibf}_{v}}$
      row_node->clust_sect_vec += row_node->V * row_node->clust_omega;
      // Expand from cluster $\cob{v}$
      for (int j = 0; j < idxs.size(); ++j) {
        y[idxs[j]] += row_node->clust_sect_vec[j];
      }
      // Recursion through row cluster tree
      ass_y_rec(row_node->sons[0]);
      ass_y_rec(row_node->sons[1]);
    }
  };
  ass_y_rec(llrcmat.rowT->root);
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
  // Not a particularly  memory-efficient implementation !
  const std::vector<HMAT::Point<1>> &row_pts{(llrcM.rowT->root)->pts};
  const std::vector<HMAT::Point<1>> &col_pts{(llrcM.colT->root)->pts};
  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < m; ++j) {
      M(i, j) = llrcM.G(row_pts[i].x[0], col_pts[j].x[0]);
    }
  }
  // Compute compressed matrix as dense matrix
  Eigen::VectorXd x(m);
  x.setZero();
  for (int j = 0; j < m; ++j) {
    x[j] = 1.0;
    if (j > 0) x[j - 1] = 0.0;
    Mt.col(j) = mvLLRPartMat(llrcM, x);
  }
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
  // Run through all far-field blocks and also determine maximal cluster count
  for (const auto &ffb : blockpart.farField) {
    // If the current block is based on a cluster, increment that cluster's
    // count in the hash map
    if (nodemap_row.find(&ffb.nx) == nodemap_row.end()) {
      nodemap_row[&ffb.nx] = 1;
      xnode_maxcnt = std::max(xnode_maxcnt, 1);
    } else {
      xnode_maxcnt = std::max(xnode_maxcnt, nodemap_row[&ffb.nx] += 1);
    }
    if (nodemap_col.find(&ffb.ny) == nodemap_col.end()) {
      nodemap_col[&ffb.ny] = 1;
      ynode_maxcnt = std::max(ynode_maxcnt, 1);
    } else {
      ynode_maxcnt = std::max(ynode_maxcnt, nodemap_col[&ffb.ny] += 1);
    }
  }
  // Run through all near-field blocks
  for (const auto &nfb : blockpart.nearField) {
    // If the current block is based on a cluster, increment that cluster's
    // count in the hash map
    if (nodemap_row.find(&nfb.nx) == nodemap_row.end()) {
      nodemap_row[&nfb.nx] = 1;
      xnode_maxcnt = std::max(xnode_maxcnt, 1);
    } else {
      xnode_maxcnt = std::max(xnode_maxcnt, nodemap_row[&nfb.nx] += 1);
    }
    if (nodemap_col.find(&nfb.ny) == nodemap_col.end()) {
      nodemap_col[&nfb.ny] = 1;
      ynode_maxcnt = std::max(ynode_maxcnt, 1);

    } else {
      ynode_maxcnt = std::max(ynode_maxcnt, nodemap_col[&nfb.ny] += 1);
    }
  }
  // Print node statistics if requested
  if (out) {
    auto output = [](const hashmap_t &nodemap, std::ostream *out) -> void {
      if (out) {
        for (const keyval_t kv : nodemap) {
          const nf_node_t *node_p = kv.first;
          const int node_count = kv.second;
          const std::vector<size_t> idxs{node_p->I()};
          (*out) << "Node with idx set = ";
          for (size_t idx : idxs) {
            (*out) << idx << ' ';
          }
          (*out) << " occurs " << node_count << " times " << std::endl;
        }
      }
    };
    (*out) << "ROW CLUSTER TREE" << std::endl;
    output(nodemap_row, out);
    (*out) << "COLUMN CLUSTER TREE" << std::endl;
    output(nodemap_col, out);
  }
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
