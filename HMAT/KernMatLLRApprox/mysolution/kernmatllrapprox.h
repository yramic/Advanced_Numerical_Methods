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
  // TO BE SUPPLEMENTED 2-4b:

  // We need to check if the input T is a cluster tree or not!
  // since T.root is already defined we know that this is an index set. 
  // Hence, it is sufficient to proof only (ii) and (iii)!

  // (ii): The subset associated with each non-leaf node is the union of the subset of its sons
  // (iii): The sets belonging to different sons of a node are disjoint

  // I want to inizialize an index Set that stores all indices of the cluster tree:
  const std::vector<std::size_t> idx_set{T.root->I()};
  // Number of points that should be contained = Max idx! Note: std::max_element returns a pointer!
  const std::size_t max_idx = *std::max_element(idx_set.begin(), idx_set.end());
  // Now I want to have a remap vec with all indices that should appear:
  std::vector<std::size_t> remap_idx(max_idx + 1, -1); // Note: it is initialized with -1 !

  for (std::size_t i{0}; i < idx_set.size(); ++i) {
    remap_idx[idx_set[i]] = i; // I basically copy the idx_set of the root!
  }

  // Recursive Checking:
  // After having initialized everything a recursive funtion needs to be created. Hence,
  // a lambda function is needed!

  //std::function returns void (only checks everything), the input is the node and it is depicted
  // as a lambda function that works by reference and takes as an input the node and the output is
  // also void!
  std::function<void(const NODE *node)> check_rec = [&] (const NODE *node) -> void {
    // First we need to check everything for the parent node!
    if (node == nullptr) {
      ok = false;
      return;
    } else {
      // Next I need a helper vector initialized with zeros which will be used later on:
      std::vector<std::size_t> helper(idx_set.size(), 0);
      // Next I want a set of indices stored from the parent tree!
      const std::vector<std::size_t> node_idx_set {node->I()};
      // This starts the check whether the indices of the parent tree are in fact contained in the root!
      for (std::size_t idx : node_idx_set) {
        if (idx > max_idx) {
          ok = false;
          return;
        } else {
          // Now I want to check if the current index is actually an index that is contained in the root!
          const std::size_t idx_check = remap_idx[idx];
          if (idx_check < 0) {
            ok = false;
            return;
          } else {
            helper[idx] += 1; // Add one on the position of the current point!
          }
        }
      }
      // Now I need to jump out of the for loop and start checking the children nodes. Since it's a binary tree
      // I need to look at two nodes and they should be disjoint from each other:
      // I have to check if the indices are contained in the parent tree and if the indices are not in the other son!
      bool has_sons = false;
      // Start Loop through sons:
      for (const NODE *son : node->sons) {
        if (son = nullptr) {
          return; // Is a leaf and has no sons!
        } else {
          has_sons = true;
          const std::vector<std::size_t> son_idx_set{son->I()};
          for (std::size_t idx_son : son_idx_set) {
            if (idx_son > max_idx) {
              // The Index is not in the root!
              ok = false;
              return;
            } else {
              const std::size_t idx_check_son = remap_idx[idx_son];
              if (idx_check_son < 0) {
                // Point not in root present!
                ok = false;
                return;
              } else {
                helper[idx_check_son] += 1;
              }
            }
          }
        }
      } // Quit for loop of sons!
      // Now I need to check the overall result whether an index appears uniquly in the child and the parent
      // or it actually appears in both children and the parent!
      if (has_sons) {
        for (std::size_t test : helper) {
          // All indices can either have the value 0 if they don't appear in a parent tree and child tree or the value
          // 2 if they appear in the parent and only ONE child not two!
          if ( (test != 0) and (test != 2) ) {
            ok = false;
            return;
          }
        } // Check over: leave for loop
      }
      // Now I actually need to start the recursion until there are no children or ok turns to false!
      if (has_sons and ok) {
        // Actually call the function for all sons!!!
        for (const NODE *son : node->sons) {
          check_rec(son);
          if (!ok) {
            // If ok becomes false return false!
            return;
          }
        }
      } 
    }
  };
 
  // Start the recursion!
  check_rec(T.root);

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
  // TO BE SUPPLEMENTED for 2-4d:

  // First I want to extract all I_row (Indices of first cluster tree)
  const std::vector<std::size_t> idx_set_row {partmat.rowT->root->I()};
  // Same for the second cluster tree:
  const std::vector<std::size_t> idx_set_col {partmat.rowT->root->I()};
  // Next, I want to determine the cardinalities of the row and column tree!
  // Find #I_x and #I_y
  const std::size_t max_row = *std::max_element(idx_set_row.begin(), idx_set_row.end());
  const std::size_t max_col = *std::max_element(idx_set_col.begin(), idx_set_col.end());
  // Now I want to initialize a helper that is used as a counter (same idea as in 2-4b!)
  // This time it's a n x m matrix:
  Eigen::Matrix<std::size_t, Eigen::Dynamic, Eigen::Dynamic> helper;
  helper.resize(idx_set_row.size(), idx_set_col.size());
  helper.setZero();
  // Next I again need a remap vector for each tree:
  std::vector<std::size_t> remap_row(idx_set_row.size(), -1);
  std::vector<std::size_t> remap_col(idx_set_col.size(), -1);
  for (std::size_t i{0}; i<remap_row.size(); ++i) {
    remap_row[idx_set_row[i]] = i;
  }
  for (std::size_t j{0}; j<remap_col.size(); ++j) {
    remap_col[idx_set_col[j]] = j;
  }
  // Example what the remap gives me:
  // I{3, 0, 2, 1} -> I_{1, 3, 2, 0} ... remap res

  // Now I need to implement two blocks of loops in the first one I go through
  // All Near Field Blocks. Whereas, in the second one the same algorithm needs
  // to be applied for Far Field Blocks!
  for (auto nfb : partmat.nearField) {
    for (std::size_t i_row : nfb.i_idx) {
      for (std::size_t i_col : nfb.j_idx) {
        if ((i_col > max_col) or (i_row > max_row)) {
          return false;
        } else {
          std::size_t check_row = remap_row[i_row];
          std::size_t check_col = remap_col[i_col];
          if ((check_col < 0) or (check_row < 0)) {
            return false;
          } else {
            helper(i_row, i_col) += 1;
          }
        }
      }
    }
  }
  // Same thing now for the Far Field Block:
  for (auto ffb : partmat.farField) {
    for (std::size_t i_row : ffb.i_idx) {
      for (std::size_t i_col : ffb.j_idx) {
        if ((i_col > max_col) or (i_row > max_row)) {
          return false;
        } else {
          std::size_t check_row = remap_row[i_row];
          std::size_t check_col = remap_col[i_col];
          if ((check_col < 0) or (check_row < 0)) {
            return false;
          } else {
            helper(i_row, i_col) += 1;
          }
        }
      }
    }
  }
  // Now I need to check all the entries of the helper, whether all of them are one
  // Note: helper.array() creates an array of the whole matrix
  bool check_helper = (helper.array() == 1).all();

  // **********************************************************************
  return check_helper;
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
<<<<<<< HEAD
// **********************************************************************
// TODO Problem 2.4e:
  // 1) After we retrieved the intervall of the boundary box [a, b], we 
  //    are now able to compute the barycentric weights Lambda from eq.
  //    eq. 2.4.6 (provided in the Hints)
  std::vector<double> lambda(q, 0); // eq. 2.4.6
  std::vector<double> t(q, 0); // Chebyshev nodes
  for (int i{1}; i < (q+1); ++i){
    const double lambda_i = (std::pow(2, q-1)/q) * std::pow(-1, i-1) * std::sin(((2*i) - 1)*M_PI/(2*q));
    lambda[i-1] = std::pow(2, q-1) * lambda_i /std::pow(b-a, q-1);
    t[i-1] = a + 0.5*(b-a) * (std::cos(((2*i) -1)*M_PI/(2*q)) + 1);
  }

  // 2) Next step is to compute the Matrix V!

  std::vector<int> helper((this->pts).size(), 1); // Stores if xi != t for the division in the end
  for (int i{0}; i < (this->pts).size(); ++i) {
    const double xi {((this->pts)[i]).x[0]}; // Collocation Point
    double tau = 0;
    for (int j{0}; j < q; ++j) {
      if (xi == t[j]) {
        helper[i] = 0;
        V(i,j) = 1; // Set to 1 for the special case according to the Lecture slide
      } else {
        tau += lambda[j] / (xi - t[j]);
        V(i,j) = lambda[j] / (xi - t[j]);
      }
    }
    if (helper[i] > 0) {
      V.row(i) /= tau;
    }
  }

// **********************************************************************
=======
    // **********************************************************************
    // TODO
    // **********************************************************************
>>>>>>> origin/master
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
<<<<<<< HEAD
// **********************************************************************
// TODO Problem 2.4f:
// First I need to define a, b and c, d the dimension of the bounding boxes
// in respect to the x and y tree:
const HMAT::BBox<1> bbox_x(_nx.pts); // Note: <1> Since the dimension needs to be 1D!
const HMAT::BBox<1> bbox_y(_ny.pts);
const double a {bbox_x.minc[0]};
const double b {bbox_x.maxc[0]};
const double c {bbox_y.minc[0]};
const double d {bbox_y.maxc[0]};
// Next, I need to compute the chebyshev nodes
std::vector<double> t_x(q), t_y(q);
for (int i{0}; i < q; ++i) {
  t_x[i] = a + 0.5*(b-a) * (std::cos((2.0*i - 1.0)*M_PI/(2.0*q)) + 1.0);
  t_y[i] = c + 0.5*(d-c) * (std::cos((2.0*i - 1.0)*M_PI/(2.0*q)) + 1.0);
}
// Next I need to fill the matrix C:
for (int k{0}; k<q; ++k) {
  for (int l{0}; l<q; ++l) {
    C(k,l) = G(t_x[k], t_y[l]);
  }
}
// **********************************************************************
=======
    // **********************************************************************
    // TODO
    // **********************************************************************
>>>>>>> origin/master
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
<<<<<<< HEAD
// Direct initialization of near field kernel collocation matrix
// **********************************************************************
// TODO Problem 2.4g:
// Since it can all be done within one for loop it is better not to store
// the collocation points and instead use them directly:
for (std::size_t i{0}; i < _nx.pts.size(); ++i) {
  for (std::size_t j{0}; j < _ny.pts.size(); ++j) {
    Mloc(i,j) = G(_nx.pts[i].x[0], _ny.pts[j].x[0]);
  }
}
// **********************************************************************
=======
  // **********************************************************************
  // TODO
  // **********************************************************************
>>>>>>> origin/master
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
// TODO Problem 2-4i:
  // The lecture slides provide a 3 pass approach where we need to use
  // the following two tools / operations:
  // (1) Restrict to Cluster
  // (2) Expand from Cluster

  // ###### PASS 1 ######
  // We want to find the vector w_w = transpose(V_w) * R_w * x
  // This has to be done as a recursion! Hence, a lambda function is needed!
  // Also this has to be done for the columne tree and not for the row tree!
  std::function<void(NODE*)> omega = [&] (NODE* n) {
    if (n == nullptr) {
      return;
    }
    // According to the hint we can use clust_omega and clust_sect_vec:
    n->clust_omega = Eigen::VectorXd::Zero(n->noIdx()); // Initialize
    n->clust_sect_vec = Eigen::VectorXd::Zero(n->noIdx()); // Initialize
    auto idx {n->I()};
    for (unsigned int i{0}; i < n->noIdx(); ++i) {
      n->clust_sect_vec(i) = x[idx[i]];
    }
    n->clust_omega = (n->V).transpose() * (n->clust_sect_vec);
    // Recursion:
    omega(n->sons[0]);
    omega(n->sons[1]);
  };
  // Start Pass 1:
  omega(llrcmat.colT->root);

  // ###### PASS 2 ######
  // First clear local storage of the Row (x) Tree!
  std::function<void(NODE*)> set_zero = [&] (NODE* n) {
    if (n == nullptr) {
      return;
    }
    n->clust_omega.setZero();
    n->clust_sect_vec.setZero();
    // Now again the recursive process down to the leafes!
    set_zero(n->sons[0]);
    set_zero(n->sons[1]);
  };
  // Start the actual recursive function process:
  set_zero(llrcmat.rowT->root);

  // Next, it is necessary to do the following two operations and differentiate
  // between far field and near field blocks:
  // 1) C x omega for the far field Blocks
  for (auto ffb : llrcmat.farField) {
    ffb.nx.clust_omega += ffb.C * ffb.ny.clust_omega;
  }
  // 2) M * R_omega for the near field Blocks:
  for (auto nfb : llrcmat.nearField) {
    nfb.nx.clust_sect_vec += nfb.Mloc * nfb.ny.clust_sect_vec;
  }

  // ###### PASS 3 ######
  // Pass 3 is now the sum of the previous passes from the rowtree:
  std::function<void(NODE*)> accumulation = [&] (NODE* n) {
    if (n == nullptr) {
      return;
    }
    n->clust_sect_vec += (n->V) * (n->clust_omega);
    auto idx = n->I();
    for (unsigned int i{0}; i < n->noIdx(); ++i) {
      y[idx[i]] += n->clust_sect_vec(i);
    }
    // Recursion down to the leafes:
    accumulation(n->sons[0]);
    accumulation(n->sons[1]);
  };
  // Start the process for pass 3:
  accumulation(llrcmat.rowT->root);
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
<<<<<<< HEAD
// **********************************************************************
  auto row_pts {(llrcM.rowT->root)->pts}; // All points from the row Tree
  auto col_pts {(llrcM.colT->root)->pts}; // All points from the column Tree
  // TODO Problem: 2-4j
  // Exact Kernel Solution:
  for (unsigned int j{0}; j < n; ++j) {
    for (unsigned int k{0}; k < m; ++k) {
      M(j,k) = llrcM.G(row_pts[j].x[0], col_pts[k].x[0]);
    }
  }
  // Solution for the compressed matrix
  // This can be achieved by a Matrix x Vector operation, look at Pass 3!
  Eigen::VectorXd x(m, 0); // Multiplication Vector initialized with 0!
  for (unsigned int i{0}; i < m; ++i) {
    x[i] = 1; // Set the current value to 1
    if (i > 0) {
      x[i-1] = 0; // Set the previous value again to zero!
    }
    Mt.col(i) = mvLLRPartMat(llrcM, x); // Matrix x Vector operation
  }
// *********************************************************************
=======
      // **********************************************************************
      // TODO
      // *********************************************************************
>>>>>>> origin/master
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
<<<<<<< HEAD
  // YOUR CODE HERE for Problem 2.4h:

  // Implementation to see if something got added!
  // std::cout << "The size of the hashmap row is: " << nodemap_row.size() << std::endl;
  // std::cout << "The size of the hashmap col is: " << nodemap_col.size() << std::endl;

  // Visualize the pairs
  // int i = 0;
  // for (const std::pair<const nf_node_t*, int>& hmat : nodemap_row) {
  //   const nf_node_t *first = hmat.first;
  //   const int counter = hmat.second;
  //   i++;
  //   std::cout << "######################################################" << std::endl;
  //   std::cout << i << "     " << *first << "     " << counter << std::endl;
  // }

  // First step: Loop over all Far field Blocks (pass by reference):
  for(auto &ffb : blockpart.farField){
    // Then, it is necessary to check whether the key exist in the hashmap. If this 
    // is not the case, this FFB have not been encountered before and we need to set
    // the value to 1! Start with the row tree:
    if (nodemap_row.find(&ffb.nx) == nodemap_row.end()) { // Note: nodemap_row.end() would be at nodemap_row.size() + 1!
      nodemap_row[&ffb.nx] = 1;
      xnode_maxcnt = std::max(xnode_maxcnt, 1); // Update if xnode_maxcnt is still 0!
    } else {
      // FFB was encountered before, hence we need to calculate plus 1 for each time it appears!
      ++nodemap_row[&ffb.nx];
      xnode_maxcnt = std::max(xnode_maxcnt, nodemap_row[&ffb.nx]); // Update xnode_maxcnt if the current value is bigger!
    }
    // Same approach for the column tree:
    if (nodemap_col.find(&ffb.ny) == nodemap_row.end()) {
      nodemap_col[&ffb.ny] = 1;
      ynode_maxcnt = std::max(ynode_maxcnt, 1);
    } else {
      ++nodemap_col[&ffb.ny];
      ynode_maxcnt = std::max(ynode_maxcnt, nodemap_col[&ffb.ny]);
    }
  }

  // Next I just copy and paste the same approach for the Near Field Block:
  for(auto &nfb : blockpart.nearField){
    // Then, it is necessary to check whether the key exist in the hashmap. If this 
    // is not the case, this FFB have not been encountered before and we need to set
    // the value to 1! Start with the row tree:
    if (nodemap_row.find(&nfb.nx) == nodemap_row.end()) { // Note: nodemap_row.end() would be at nodemap_row.size() + 1!
      nodemap_row[&nfb.nx] = 1;
      xnode_maxcnt = std::max(xnode_maxcnt, 1); // Update if xnode_maxcnt is still 0!
    } else {
      // FFB was encountered before, hence we need to calculate plus 1 for each time it appears!
      ++nodemap_row[&nfb.nx];
      xnode_maxcnt = std::max(xnode_maxcnt, nodemap_row[&nfb.nx]); // Update xnode_maxcnt if the current value is bigger!
    }
    // Same approach for the column tree:
    if (nodemap_col.find(&nfb.ny) == nodemap_row.end()) {
      nodemap_col[&nfb.ny] = 1;
      ynode_maxcnt = std::max(ynode_maxcnt, 1);
    } else {
      ++nodemap_col[&nfb.ny];
      ynode_maxcnt = std::max(ynode_maxcnt, nodemap_col[&nfb.ny]);
    }
  }


  // Through looking at the solution I just saw that for our case it would be sufficient to only compute:
  // ++nodemap_row[&ffb.nx] and ++nodemap_col[&ffb.ny] same for NFB! This would save at least some code, but
  // I am not sure if this approach would be in general correct! Hence I sticked to this approach!

  // Also it is possible not to update xnode_maxcnt and ynode_maxcnt at each step, it would be sufficient to append
  // the following line in the end (Note: Same approach for ynode_maxcnt):
  // xnode_maxcnt = std::max_element(nodemap_row.begin(), nodemap_row.end(),
  //                                 key_less)->second;


  // The following code checks the result!
  // int i = 0;
  // for (const std::pair<const nf_node_t*, int>& hmat : nodemap_row) {
  //   const nf_node_t *first = hmat.first;
  //   const int counter = hmat.second;
  //   i++;
  //   std::cout << "######################################################################" << std::endl;
  //   std::cout << i << "     " << *first << "     " << counter << std::endl;
  // }

  // std::cout << "The size of the hashmap row is: " << nodemap_row.size() << std::endl;
  // std::cout << "The size of the hashmap col is: " << nodemap_col.size() << std::endl;

// **********************************************************************
=======
  // YOUR CODE HERE
  // **********************************************************************
>>>>>>> origin/master
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
