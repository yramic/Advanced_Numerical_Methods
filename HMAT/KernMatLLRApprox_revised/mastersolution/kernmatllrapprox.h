/**
 * @file kernmatllrapprox.h
 * @brief NPDE homework KernMatLLRApprox code
 * @author R. Hiptmair, Peiyuan Xie
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
  const std::vector<std::size_t> idxset{T.root->I};
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
    const std::vector<std::size_t> idxset_node{node->I};  // Index set of node
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
    for (auto &son : node->sons) {
      if (son != nullptr) {
        has_sons = true;
        const std::vector<std::size_t> idxset_son{son->I}; // son's index set 
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
        for (auto &son : node->sons) {
	  ctcheck_rec(son.get());
	  if (!ok) return;
	} } }
  };
  // Lanch recursion
  ctcheck_rec(T.root.get());
  return ok;
}
/* SAM_LISTING_END_1 */
// clang-format on

// clang-format off
/** @brief Verify validity of matrix paving into near and far field blocks */
/* SAM_LISTING_BEGIN_2 */
template <typename NODE>
bool checkMatrixPartition(
    const HMAT::BlockPartition<NODE> &partmat) {
  assertm((partmat.rowT and partmat.colT), "Missing trees!");
  // Index sets $\mathbb{I}$, $\mathbb{J}$ underlying the row/col cluster trees
  const std::vector<std::size_t> idxset_row{partmat.rowT->root->I};
  const std::vector<std::size_t> idxset_col{partmat.colT->root->I};
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

/** Extended class for cluster trees for local low-rank approximation */
/* SAM_LISTING_BEGIN_E */
class LLRClusterTree : public HMAT::ClusterTree<HMAT::CtNode<1>> {
 public:
  // Idle constructor just setting rank argument q
  explicit LLRClusterTree(size_t _q, const std::vector<HMAT::Point<1>> &pts,
                          std::size_t minpts = 1)
      : HMAT::ClusterTree<HMAT::CtNode<1>>(pts, minpts), q(_q) {
    // Resize global vectors
    Vs.resize(this->ptsT.size(), q);
    clust_omega.resize(q, this->numNodes);
    clust_sect_vec.resize(this->ptsT.size());

    // Recursive initialization of the low-rank factor matrices \cob{$\VV_w$}
    std::function<void(HMAT::CtNode<1> * node)> initvrec =
        [&](HMAT::CtNode<1> *node) -> void {
      // Retrieve bounding box,
      // \href{https://stackoverflow.com/questions/4643074/why-do-i-have-to-access-template-base-class-members-through-the-this-pointer}{explanation}
      // for the this pointer
      const HMAT::BBox<1> bbox = this->getBBox(node);
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
        lambda[i] = fac * sgn * std::sin(arg);  // \prbeqref{eq:bwf}
      }
      // Number of collocation points in cluster
      const unsigned int nIw = node->noIdx();
      // Traverse collocation points (in the interval $\cintv{a,b}$)
      Eigen::VectorXd tx_diff(q);  // $t_i - \xi$
      auto V = getV(node);
      for (int j = 0; j < nIw; ++j) {
        // Coordinate of current collocation point
        const double xi = ((this->ptsT)[node->offset + j]).x[0];
        double tau =
            0.0;  // $\cob{\tau := \sum_{i=1}^{q} \frac{\lambda_i}{\xi-t_i}}$
        bool on_node = false;
        for (int i = 0; i < q; i++) {
          tx_diff[i] = xi - t[i];
          // Avoid division by zero
          if (tx_diff[i] == 0.0) {
            on_node = true;
            V.row(j).setZero();
            V(j, i) = 1.0;  // $\cob{(\VV)_{j,:} = \Ve_{i}^{\top}}$ when hitting
                            // a node
            break;          // The $j$-th row of $\VV$ is complete already
          } else {
            const double txdl = lambda[i] / tx_diff[i];
            V(j, i) = txdl;
            tau += txdl;
          }
        }
        if (!on_node) V.row(j) /= tau;
      }

      if (node->sons[0]) initvrec(node->sons[0].get());
      if (node->sons[1]) initvrec(node->sons[1].get());
    };
    initvrec(this->root.get());
  }
  virtual ~LLRClusterTree() = default;

 public:
  // For the conciseness of access
  Eigen::VectorBlock<Eigen::VectorXd> getSectVec(const HMAT::CtNode<1> *node) {
    return clust_sect_vec.segment(node->offset, node->noIdx());
  }
  Eigen::MatrixXd::ColXpr getOmega(const HMAT::CtNode<1> *node) {
    return clust_omega.col(node->nodeNumber);
  }
  Eigen::Block<Eigen::MatrixXd> getV(const HMAT::CtNode<1> *node) {
    return Vs.middleRows(node->offset, node->noIdx());
  }

  // rank of separable approximation on cluster boxes big matrix
  const std::size_t q;
  // global matrix of low-rank factors, each \cob{$\VV\in\bbR^{k,q}$}
  Eigen::MatrixXd Vs;
  // global matrix for cluster-local linear algebra
  Eigen::MatrixXd clust_omega;
  // temporary storage for cluster-associated vector sections
  Eigen::VectorXd clust_sect_vec;
};
/* SAM_LISTING_END_E */

/** Extended class for block partition, knowing low-rank compression */
/* SAM_LISTING_BEGIN_H */
template <class TREE, typename KERNEL>
class BiDirChebBlockPartition : public HMAT::BlockPartition<TREE> {
 public:
  BiDirChebBlockPartition(std::shared_ptr<TREE> _rowT,
                          std::shared_ptr<TREE> _colT, KERNEL _Gfun,
                          std::size_t _q, double eta0 = 2.0)
      : HMAT::BlockPartition<TREE>(_rowT, _colT), G(_Gfun), q(_q) {
    HMAT::BlockPartition<TREE>::init(eta0);

    // build C matrix for every far field block
    Cs.resize(this->farField.size());
    for (int n = 0; n < Cs.size(); n++) {
      // get nodes in the block
      NODE *nx = &this->farField[n].nx;
      NODE *ny = &this->farField[n].ny;
      // resize the matrix for this block
      Cs[n].resize(q, q);
      // Obtain bounding boxes, here intervals
      std::array<HMAT::BBox<1>, 2> bboxs = {this->rowT->getBBox(nx),
                                            this->colT->getBBox(ny)};
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
          t(d, j) =
              intv[d][0] + 0.5 * (intv[d][1] - intv[d][0]) * (cosval + 1.0);
        }
      }
      // Fill matrix $\cob{\VC}$: $\cob{(\VC)_{k,\ell} = G(t_k,t_{\ell})}$, see
      // \prbeqref{eq:Cm}
      for (int k = 0; k < q; ++k) {
        for (int j = 0; j < q; ++j) {
          Cs[n](k, j) = G(t(0, k), t(1, j));
        }
      }
    }

    // build M matrix for every near field block
    Mlocs.resize(this->nearField.size());
    for (int n = 0; n < Mlocs.size(); n++) {
      // get nodes in this block
      NODE *nx = &this->nearField[n].nx;
      NODE *ny = &this->nearField[n].ny;
      // resize the matrix of this block
      Mlocs[n].resize(nx->noIdx(), ny->noIdx());
      // Direct initialization of near field kernel collocation matrix
      for (int i = 0; i < Mlocs[n].rows(); ++i) {
        for (int j = 0; j < Mlocs[n].cols(); ++j) {
          Mlocs[n](i, j) = G((this->rowT->ptsT[nx->offset + i]).x[0],
                             (this->colT->ptsT[ny->offset + j]).x[0]);
        }
      }
    }
  }
  virtual ~BiDirChebBlockPartition() = default;

 protected:
  using NODE = typename TREE::node_t;
  // Construct an instance of far-field block type
  virtual HMAT::IndexBlock<NODE> makeFarFieldBlock(NODE &nx, NODE &ny) {
    ffb_cnt++;
    return HMAT::IndexBlock<NODE>(nx, ny);
  }
  // Construct an instance of near-field block type
  virtual HMAT::IndexBlock<NODE> makeNearFieldBlock(NODE &nx, NODE &ny) {
    nfb_cnt++;
    return HMAT::IndexBlock<NODE>(nx, ny);
  }

 public:
  KERNEL G;                 // Reference to the kernel function
  const std::size_t q;      // degree+1 of interpolating polynomial
  unsigned int ffb_cnt{0};  // Counter for far-field blocks
  unsigned int nfb_cnt{0};  // Counter for near-field blocks
  std::vector<Eigen::MatrixXd>
      Cs;  // vector containing matrices C for every far field block
  std::vector<Eigen::MatrixXd>
      Mlocs;  // vector containing matrices M for every near field block
};
/* SAM_LISTING_END_H */

// Special data type for local low-rank compression by one-dimensional
// bi-directional Chebychev interpolation
template <typename KERNEL>
using BiDirChebPartMat1D = BiDirChebBlockPartition<LLRClusterTree, KERNEL>;

// clang-format off
// Matrix x Vector based on compressed kernel collocation matrix
// Assumes contiguous indices in root node!
// The NODE template parameter must be compatible with InterpNode
/* SAM_LISTING_BEGIN_U */
template <class TREE, typename KERNEL>
Eigen::VectorXd mvLLRPartMat(BiDirChebBlockPartition<TREE, KERNEL> &llrcmat,
                             const Eigen::VectorXd &x) {
  using NODE = typename TREE::node_t;
  // Verify requirement of contiguous indices
  const std::vector<size_t> col_idxs = (llrcmat.colT->root)->I;
  assertm((*std::max_element(col_idxs.begin(), col_idxs.end()) ==
           col_idxs.size() - 1),
          "col idxs not contiguous");
  const std::vector<size_t> row_idxs = (llrcmat.rowT->root)->I;
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
      const std::vector<size_t> idxs = col_node->I;
      // Restriction of argument vector to column cluster, assuming contiguous indices in idxs
      llrcmat.colT->getSectVec(col_node) = x.segment(idxs.front(), idxs.size());
      // Compute $\cob{\vec{\omegabf}_{w}}$
      auto V = llrcmat.colT->getV(col_node);
      Eigen::VectorXd Rmu = llrcmat.colT->getSectVec(col_node);
      llrcmat.colT->getOmega(col_node) = V.transpose() * Rmu;
      // Recursion: visit entire cluster tree
      comp_omega_rec(col_node->sons[0].get());
      comp_omega_rec(col_node->sons[1].get());
    }
  };
  comp_omega_rec(llrcmat.colT->root.get());
  // Final sentence in \lref{par:3p}: Clear local storage of row tree
  std::function<void(NODE * row_node)> clear_vec_rec = [&](NODE *row_node) -> void {
    if (row_node) {
      llrcmat.rowT->getSectVec(row_node).setZero();
      llrcmat.rowT->getOmega(row_node).setZero();
      clear_vec_rec(row_node->sons[0].get());
      clear_vec_rec(row_node->sons[1].get());
    }
  };
  clear_vec_rec(llrcmat.rowT->root.get());

  // Pass (II) in \lref{par:3p}: Do block-based computations
  // Visit far-field blocks:
  for (int n=0; n < llrcmat.farField.size(); n++) {
    // Far-field block: accumulate $\cob{\VC_{v\times w}\vec{\omegabf}_w}$
    NODE &row_node = llrcmat.farField[n].nx;
    const NODE &col_node = llrcmat.farField[n].ny;
    llrcmat.rowT->getOmega(&row_node) += llrcmat.Cs[n] * llrcmat.colT->getOmega(&col_node);
  }
  for (int n=0; n < llrcmat.nearField.size(); n++) {
    // Near-field block: use retricted kernel collocatin matrix $\cob{\rst{\VM}{v\times w}}$
    NODE &row_node = llrcmat.nearField[n].nx;
    const NODE &col_node = llrcmat.nearField[n].ny;
    llrcmat.rowT->getSectVec(&row_node) += llrcmat.Mlocs[n] * llrcmat.colT->getSectVec(&col_node);
  }

  // Pass (III) of \lref{par:3p}: Traverse row tree and assemble return vector
  std::function<void(NODE * row_node)> ass_y_rec = [&](NODE *row_node) -> void {
    if (row_node) {
      const std::vector<size_t> idxs = row_node->I;
      // $\cob{\rst{\Vx}{v} += \VU_{v}\vec{\zetabf}_{v}+\vec{\phibf}_{v}}$
      auto V = llrcmat.rowT->getV(row_node);
      auto omega = llrcmat.rowT->getOmega(row_node);
      llrcmat.rowT->getSectVec(row_node) += V * omega;
      // Expand from cluster $\cob{v}$, assuming contiguous indices in idxs
      y.segment(idxs.front(), idxs.size()) += llrcmat.rowT->getSectVec(row_node);
      // Recursion through row cluster tree
      ass_y_rec(row_node->sons[0].get());
      ass_y_rec(row_node->sons[1].get());
    }
  };
  ass_y_rec(llrcmat.rowT->root.get());
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
  const std::vector<HMAT::Point<1>> &row_pts{(llrcM.rowT)->ptsT};
  const std::vector<HMAT::Point<1>> &col_pts{(llrcM.colT)->ptsT};
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
template <typename TREE>
unsigned int computeSparsityMeasure(const HMAT::BlockPartition<TREE> &blockpart,
                                    std::ostream *out = nullptr) {
  assertm((blockpart.rowT and blockpart.colT), "Missing trees!");
  // Set up hash maps for nodes of both trees
  using NODE = typename TREE::node_t;
  using hashmap_t = std::unordered_map<const NODE *, int>;
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
          const NODE *node_p = kv.first;
          const int node_count = kv.second;
          const std::vector<size_t> idxs{node_p->I};
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
