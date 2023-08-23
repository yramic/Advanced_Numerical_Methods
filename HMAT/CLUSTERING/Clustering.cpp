/***********************************************************************
 *                                                                     *
 * Code for Course "Advanced Numerical Methods for CSE"                *
 * (Prof. Dr. R. Hiptmair)                                             *
 * Author: R.H.                                                        *
 * Date: Nov 18, 2017                                                  *
 * (C) Seminar for Applied Mathematics, ETH Zurich                     *
 * This code can be freely used for non-commercial purposes as long    *
 * as this header is left intact.                                      *
 ***********************************************************************/

// General includes
#include <Eigen/Dense>
#include <algorithm>
#include <cmath>
#include <exception>
#include <iostream>
#include <limits>
#include <utility>
#include <vector>

namespace HMAT {


// ======================================================================
// For local low-rank approximation
// ======================================================================

/** Node supporting separable approximation by interpolation */
/* SAM_LISTING_BEGIN_Y */
template <int d>
class InterpNode : public CtNode<d> {
 public:
  // Constructor from sequence of points; initializes \cob{$\VV$}
  InterpNode(const std::vector<Point<d>> _pts, std::size_t _q, int _dir = 0)
      : CtNode<d>(_pts, _dir), q(_q), sons{nullptr, nullptr}, k(_pts.size()) {
    initV();
  }
  virtual ~InterpNode() {}

 protected:
  // Initialization of matrix \cob{$\VV$}
  void initV();

 public:
  const int q;                       // Rank, no of interpolation nodes
  std::size_t k;                     // Number of indices contained
  Eigen::MatrixXd V;                 // low-rank factor \cob{$\VV\in\bbR^{k,q}$}
  std::array<InterpNode *, 2> sons;  // Pointers to sons (of type InterpNode!)
};
/* SAM_LISTING_END_Y */

// >> Implementation required

/* SAM_LISTING_BEGIN_Z */
template <int d>
void InterpNode<d>::initV() {
  static_assert(d == 1, "Implemented only for 1D");
  std::cerr << "Rank-" << q
            << " InterpNode: Initialization of V not implemented" << std::endl;
}

/* SAM_LISTING_END_Z */

template <int dim>
std::ostream &operator<<(std::ostream &o, const InterpNode<dim> &node) {
  o << "## IPNode: " << node.pts.size() << " points, " << node.getBBox()
    << ": ";
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
class LowRankClusterTree : public ClusterTree<Node> {
 public:
  // Idle constructor just setting rank argument q
  explicit LowRankClusterTree(size_t _q) : q(_q) {}
  // Actual constructor taking a sequence of points
  void init(const std::vector<Point<Node::dim>> pts, std::size_t minpts = 1);
  virtual ~LowRankClusterTree() {}

 protected:
  // factory method for relevant type of node taking rank argument
  virtual Node *createNode(const std::vector<Point<Node::dim>> pts, int dir) {
    return new Node(pts, q, dir);
  }

 public:
  const std::size_t q;  // rank of separable approximation on cluster boxes
};

template <class Node>
void LowRankClusterTree<Node>::init(const std::vector<Point<Node::dim>> pts,
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
  BiDirChebInterpBlock(const Node &nx, const Node &ny, kernel_t _G,
                       std::size_t _q);
  virtual ~BiDirChebInterpBlock() {}
  // Invalid constructor throwing exception
  BiDirChebInterpBlock(const Node &nx, const Node &ny);

  const kernel_t G;   // kernel function \cob{$\krn$}
  const int q;        // No of interpolation nodes
  Eigen::MatrixXd C;  // \cob{$\VC\in\bbR^{q,q}$}
};
/* SAM_LISTING_END_F */

template <class Node, typename KERNEL>
BiDirChebInterpBlock<Node, KERNEL>::BiDirChebInterpBlock(const Node &_nx,
                                                         const Node &_ny,
                                                         KERNEL _G,
                                                         std::size_t _q)
    : IndexBlock<Node>(_nx, _ny), G(_G), q(_q) {
  std::cerr << "BiDirChebInterpBlock: q = " << q
            << ": Initialization not yet implemented" << std::endl;
}

template <class Node, typename KERNEL>
BiDirChebInterpBlock<Node, KERNEL>::BiDirChebInterpBlock(const Node &_nx,
                                                         const Node &_ny)
    : IndexBlock<Node>(_nx, _ny), G(kernel_t()), q(0) {
  throw std::runtime_error(
      "BiDirChebInterpBlock: no construction without kernel");
}

/** General type for generic near-field cluster pair */
/* SAM_LISTING_BEGIN_G */
template <class Node, typename KERNEL>
class NearFieldBlock : public IndexBlock<Node> {
 public:
  using kernel_t = KERNEL;
  NearFieldBlock(const Node &nx, const Node &ny, kernel_t _G);
  virtual ~NearFieldBlock() {}
  // Invalid constructor throwing exception
  NearFieldBlock(const Node &nx, const Node &ny);

  const kernel_t G;      // kernel function \cob{$\krn$}
  Eigen::MatrixXd Mloc;  // local kernel collocation matrix
};
/* SAM_LISTING_END_G */

template <class Node, typename KERNEL>
NearFieldBlock<Node, KERNEL>::NearFieldBlock(const Node &_nx, const Node &_ny,
                                             kernel_t _G)
    : IndexBlock<Node>(_nx, _ny), G(_G) {
  std::cerr << "NearFieldBlock: Initialization not yet implemented"
            << std::endl;
}

template <class Node, typename KERNEL>
NearFieldBlock<Node, KERNEL>::NearFieldBlock(const Node &_nx, const Node &_ny)
    : IndexBlock<Node>(_nx, _ny), G(kernel_t()) {
  throw std::runtime_error("Near Field Block: no construction without kernel");
}

/** Extended class for block partition, knowing low-rank compression */
/* SAM_LISTING_BEGIN_H */
template <class Node, typename FFB, typename NFB>
class HierMatBlockPartition : public BlockPartition<Node, FFB, NFB> {
 public:
  using kernel_t = typename NFB::kernel_t;
  HierMatBlockPartition(const LowRankClusterTree<Node> &_xT,
                        const LowRankClusterTree<Node> &_yT, kernel_t _G,
                        std::size_t _q, double eta0 = 0.5)
      : BlockPartition<Node, FFB, NFB>(_xT, _yT), G(_G), q(_q) {
    BlockPartition<Node, FFB, NFB>::init(eta0);
  }
  virtual ~HierMatBlockPartition() {}

 protected:
  // Construct an instance of far-field block type
  virtual FFB makeFarFieldBlock(const Node &nx, const Node &ny) const {
    std::cout << "HierMatBlockPartition: makeFarFieldBlock" << std::endl;
    return FFB(nx, ny, G, q);
  }
  // Construct an instance of near-field block type
  virtual NFB makeNearFieldBlock(const Node &nx, const Node &ny) const {
    std::cout << "HierMatBlockPartition: makeNearFieldBlock" << std::endl;
    return NFB(nx, ny, G);
  }

 public:
  KERNEL G;     // Reference to the kernel function
  const std::size_t q;  // degree+1 of interpolating polynomial
};
/* SAM_LISTING_END_H */

/** Class providing the kernel function */
struct Kernel {
  double operator()(double x, double y) {
    if (x != y)
      return -std::log(std::abs(x - y));
    else
      return 0.0;
  }
};

/* SAM_LISTING_BEGIN_R */
template <int d>  // dimension \cob{$d$} as template argument
struct BasisFn {
  std::size_t idx;                         // Index of basis function
  Eigen::Matrix<double, d, 1> xmin, xmax;  // Corners of bounding box
};
/* SAM_LISTING_END_R */

}  // end namespace HMAT

/** Main program */
int main(int, char **) {
  const int d = 1;
  std::size_t npts = 16;
  std::vector<HMAT::Point<d>> pts;
  for (int n = 0; n < npts; n++) {
    HMAT::Point<d> p;
    p.idx = n;
    p.x[0] = (double)n;
    pts.push_back(p);
  }
  {
    HMAT::ClusterTree<HMAT::CtNode<d>> T;
    T.init(pts);
    std::cout << T << std::endl;

    HMAT::BlockPartition<HMAT::CtNode<d>, HMAT::IndexBlock<HMAT::CtNode<d>>,
                         HMAT::IndexBlock<HMAT::CtNode<d>>>
        bP(T, T);
    bP.init();
    std::cout << bP << "Part I: Normal termination" << std::endl;
  }
  {
    std::size_t q = 7;
    HMAT::Kernel G;
    HMAT::LowRankClusterTree<HMAT::InterpNode<d>> LRT(q);
    LRT.init(pts);
    std::cout << "LowRankClusterTree:" << std::endl << LRT << std::endl;
    HMAT::HierMatBlockPartition<
        HMAT::InterpNode<d>,
        HMAT::BiDirChebInterpBlock<HMAT::InterpNode<d>, HMAT::Kernel>,
        HMAT::NearFieldBlock<HMAT::InterpNode<d>, HMAT::Kernel>>
        lrbP(LRT, LRT, G, q);
    std::cout << "Part II: Normal termination" << std::endl;
  }
  exit(0);
}
