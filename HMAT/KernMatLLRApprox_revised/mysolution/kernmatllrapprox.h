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
    const HMAT::BlockPartition<NODE> &partmat) {
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

// **********************************************************************
// TODO
// **********************************************************************
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
  // vector containing matrices C for every far field block
  std::vector<Eigen::MatrixXd> Cs;
  // vector containing matrices M for every near field block
  std::vector<Eigen::MatrixXd> Mlocs;
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
