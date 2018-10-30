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
#include <iostream>
#include <cmath>
#include <limits>
#include <vector>
#include <utility>
#include <algorithm>
#include <exception>
#include <Eigen/Dense>

using namespace std;
using namespace Eigen;

/** @brief Data structure for a collocation point
   A collocation point has an index and coordinates 
   @param t dimension of ambient space */
/* SAM_LISTING_BEGIN_1 */
template<int d>  // dimension \cob{$d$} as template argument
struct Point {
  size_t idx;            // number of collocation point
  Matrix<double,d,1> x;  // coordinate vector
};
/* SAM_LISTING_END_1 */

/** @brief Data structure for bounding box
   a bounding box is defined by two position vectors */
/* SAM_LISTING_BEGIN_2 */
template<int d>  // dimension \cob{$d$} as template argument
struct BBox {
  // Bounding box from sequence of points
  BBox(const vector<Point<d>> pts);
  // Size \cob{$\diam(B)$} of a bounding box
  double diam(void) const {
    return (maxc-minc).cwiseAbs().maxCoeff(); }
  // Corner points of bounding box
  Matrix<double,d,1> minc,maxc;
};
// distance of \cob{$\cintv{a,b}$} and \cob{$\cintv{c,d}$}
double dist(double a,double b,double c,double d) {
  if (b<a) swap(a,b); if (d<c) swap(c,d);
  if (c<a) { swap(a,c); swap(b,d); }
  return (c<b)?0.0:c-b;
}
// distance of d-dimensional boxes
template<int d>
double dist(const BBox<d> &bx,const BBox<d> &by) {
  double dst = 0.0;
  for(int l=0;l<d;++l) 
    dst += pow(dist(bx.minc[l],bx.maxc[l],by.minc[l],by.maxc[l]),2);
  return sqrt(dst);
}
/* SAM_LISTING_END_2 */

//
template<int d>
ostream &operator << (ostream &o,const BBox<d> &box) {
  return o << "BBOX: " << box.minc << ',' << box.maxc << ' ';
}

/* SAM_LISTING_BEGIN_3 */
template<int d>
BBox<d>::BBox(const vector<Point<d>> pts) {
  double tmp;
  minc = numeric_limits<double>::max()*Matrix<double,d,1>::Ones();
  maxc = -minc;
  for(const Point<d> &v: pts) {
    for(int l=0;l<d;l++) {
      const double tmp = v.x[l];
      if (tmp < minc[l]) minc[l] = tmp;
      if (tmp > maxc[l]) maxc[l] = tmp;
    }}
} // end BBox constructor
/* SAM_LISTING_END_3 */

/** @brief Data structure for the node of a binary cluster tree
   A node of a cluster tree contains a set of collocation points. */
/* SAM_LISTING_BEGIN_4 */
template <int d>
class CtNode {
public:
  const static size_t dim = d;
  // Constructors taking a sequence of points
  CtNode(const vector<Point<d>> _pts,int _dir=0):
    pts(_pts),sons{nullptr,nullptr},dir(_dir) {}
  // Destructor (also attempts to destroy sons!)
  virtual ~CtNode(void) {
    if (sons[0]) delete sons[0];
    if (sons[1]) delete sons[1];
  }
  // Number of indices owned by the cluster
  size_t noIdx(void) const { return pts.size(); } 
  // Function \cob{$\Ci$}: access to owned indices
  vector<size_t> I(void) const;
  // Access to bounding box
  BBox<d> getBBox(void) const { return BBox<d>(pts); }
  // Is the node a leaf node ?
  inline bool isLeaf(void) const { return(!sons[0] || !sons[1]); }
  // Output operator
  template<int dim>
  friend  ostream &operator << (ostream &o,const CtNode<dim> &node); 
  // Pointers to sons
  CtNode *sons[2];
  // Points contained in the cluster
  vector<Point<d>> pts;
  // Direction for sorting
  int dir;
};
/* SAM_LISTING_END_4 */


/* SAM_LISTING_BEGIN_8 */
template<int d>
vector<size_t> CtNode<d>::I(void) const {
  vector<size_t> idx;
  for (const Point<d> &pt: pts) idx.push_back(pt.idx);
  return idx;
}
/* SAM_LISTING_END_8 */

template<int dim>
ostream &operator << (ostream &o,const CtNode<dim> &node) {
  o << "Node: " << node.pts.size()  << " points, " << node.getBBox() << ": ";
  for(const Point<dim> &pt: node.pts) {
    o << "["; for(int l=0;l<dim;l++) { o << pt.x[l] << ' '; } o << "]";
  }
  o << endl;
  if (node.sons[0]) o << *(node.sons[0]);
  if (node.sons[1]) o << *(node.sons[1]);
  return o;
}

/** @brief Data structure for a cluster tree.
    A cluster tree builds and describes a multilevel partitioning of a set of points */
/* SAM_LISTING_BEGIN_5 */
template <class Node>
class ClusterTree {
public:
  const static size_t dim = Node::dim; // space dimension d
  // Idle constructor
  ClusterTree(void):root(nullptr) {}
  // Effective Constructor taking a sequence of points
  // (needed, because polynorphism not supported in constructor)
  void init(const vector<Point<dim>> pts,size_t minpts = 1);
  virtual ~ClusterTree(void) { if (root) delete root; }
  // Output of tree
  template <class Nd>
  friend ostream &operator << (ostream &o,const ClusterTree<Nd> &T);
protected:
  // Recursive construction
  virtual void buildRec(Node *nptr,size_t minpts);
  // Node factory
  virtual Node* createNode(const vector<Point<dim>> pts,int dir) {
    return new Node(pts,dir); }
 public:
  Node *root; // pointer to root node
};
/* SAM_LISTING_END_5 */

/* SAM_LISTING_BEGIN_6 */
template<class Node>
void ClusterTree<Node>::init(const vector<Point<dim>> pts,size_t minpts) {
  if (!(root = createNode(pts,0)))
    throw(runtime_error("Cannot allocate root"));
  if (minpts < 1)
    throw(runtime_error("minpts must be at least 1"));
  buildRec(root,minpts);
}
/* SAM_LISTING_END_6 */

template <class Node>
ostream &operator << (ostream &o,const ClusterTree<Node> &T) {
  return o << "ROOT: " << *(T.root) << "END CLUSTER TREE" << endl;
}
  
/* SAM_LISTING_BEGIN_7 */
template<class Node>
void ClusterTree<Node>::buildRec(Node *nptr,size_t minpts) {
  const size_t n = nptr->pts.size(); // Number of held indices
  // Leaf, if minimal number of indices reached
  if (n > minpts) { // \Label[line]{brc:1}
    // Points have to be copied and sorted according to direction dir
    vector<Point<dim>> tpts(nptr->pts);
    // next sorting direction 
    const int dir = (nptr->dir + 1)%dim;
    // call sort function from standard library
    sort(tpts.begin(),tpts.end(),
	 [dir] (const Point<dim> &p1,const Point<dim> &p2)
	 -> bool { return (bool)(p1.x[dir] < p2.x[dir]); });
    // Split point sequence and construct sons
    const size_t m = n/2; // integer arithmeric, m>0 ensured
    const vector<Point<dim>> low_pts(tpts.cbegin(),tpts.cbegin()+m);
    // First son gets ``lower half'' of sorted points
    if (!(nptr->sons[0] = createNode(low_pts,dir)))
      throw(runtime_error("Cannot allocate first son"));
    buildRec(nptr->sons[0],minpts); // recurse into first son
    // Second son get ``upper half'' of sorted points
    const vector<Point<dim>> up_pts(tpts.cbegin()+m,tpts.cend());
    if (!(nptr->sons[1] = createNode(up_pts,dir)))
      throw(runtime_error("Cannot allocate second son"));
    buildRec(nptr->sons[1],minpts); // recurse into 2nd son
  }
}
/* SAM_LISTING_END_7 */

/** @brief Data structure for both far-field and near-field blocks */
/* SAM_LISTING_BEGIN_9 */
template <class Node>
struct IndexBlock {
  // Constructors extracts indices from clusters
  IndexBlock(const Node &_nx,const Node &_ny):
    nx(_nx),ny(_ny),i_idx(_nx.I()),j_idx(ny.I()) {}
  virtual ~IndexBlock(void) {}
  const Node &nx,&ny;               // contributing clusters
  const vector<size_t> i_idx,j_idx; // contained indices
};
/* SAM_LISTING_END_9 */

/* SAM_LISTING_BEGIN_A */
template <class Node,typename FFB,typename NFB>
class BlockPartition {
public:
  // Idle constructor
  BlockPartition(const ClusterTree<Node> &_xT,
		 const ClusterTree<Node> &_yT):
    xT(_xT),yT(_yT) {}
  // Trigger recursive construction of partition
  // (Needed, because polymorphic functions not available in constructor)
  void init(double eta0 = 0.5); 
  virtual ~BlockPartition(void) {}
  // Admissibility condition \cob{$\adm$}, see \cref{def:ac}
  virtual bool adm(const Node *nx,const Node *ny,double eta0) const;
protected:
  // Recursive construction from cluster pair
  virtual void buildRec(const Node *nx,const Node *ny,double eta0);
  // Construct an instance of far-field block type
  virtual FFB makeFarFieldBlock(const Node &nx,const Node &ny) const
  { return FFB(nx,ny); }
  // Construct an instance of near-field block type
  virtual NFB makeNearFieldBlock(const Node &nx,const Node &ny) const
  { return NFB(nx,ny); }
public:
  const ClusterTree<Node> &xT,&yT; // underlying cluster trees
  vector<FFB> farField;  // index blocks in the far field
  vector<NFB> nearField; // index blocks in the near field
  static bool dbg; // Debugging flag
};
/* SAM_LISTING_END_A */

template <class Node,typename FFB,typename NFB>
bool BlockPartition<Node,FFB,NFB>::dbg = false;

/* SAM_LISTING_BEGIN_B */
template <class Node,typename FFB,typename NFB>
void BlockPartition<Node,FFB,NFB>::init(double eta0) {
  buildRec(xT.root,yT.root,eta0);
}
/* SAM_LISTING_END_B */

/* SAM_LISTING_BEGIN_C */
template <class Node,typename FFB,typename NFB>
bool BlockPartition<Node,FFB,NFB>::adm(const Node *nx,
				       const Node *ny,
				       double eta0) const {
  // Neither node must be a leaf.
  if (nx->isLeaf() || ny->isLeaf()) return false;
  // Geometric admissibility condition, see \cref{eq:etadef}.
  const BBox<Node::dim> Bx = nx->getBBox(),By = ny->getBBox();
  const double eta = max(Bx.diam(),By.diam())/(2*dist(Bx,By));
  return (eta < eta0);
}
/* SAM_LISTING_END_C */

/* SAM_LISTING_BEGIN_D */
template <class Node,typename FFB,typename NFB>
void BlockPartition<Node,FFB,NFB>::
buildRec(const Node *nx,const Node *ny,double eta0) {
  if (nx && ny) {
    // Add admissible pair to far field
    if (adm(nx,ny,eta0)) // \Label[line]{bpr:adm}
      farField.push_back(makeFarFieldBlock(*nx,*ny));
    else {
      bool rec = false;
      for (int isx=0;isx <=1; isx++) { 
	for (int isy=0;isy <=1; isy++) {
	  if (nx->sons[isx] && ny->sons[isy]) {
	    // Next level of recursion for non-leaves
	    rec= true; buildRec(nx->sons[isx],ny->sons[isy],eta0);
	  }}}
      // At least one leaf cluster:
      // Add cluster pair to near field
      if (!rec) // \Label[line]{bpr:1}
	nearField.push_back(makeNearFieldBlock(*nx,*ny));
    }
  }
  else 
    throw(runtime_error("Invalid node pointers"));
}
/* SAM_LISTING_END_D */

template <class Node,typename FFB,typename NFB>
ostream &operator << (ostream &o,const BlockPartition<Node,FFB,NFB> &bp) {
  o << "Near field indices:" << endl;
  for(const NFB &b: bp.nearField) {
    o << "[ i_idx = "; for(int i: b.i_idx) o << i << ", "; o << endl;
    o << "  j_idx = "; for(int j: b.j_idx) o << j << ", "; o << "]" << endl;
  }
  o << "Far field indices:" << endl;
  for(const FFB &b: bp.farField) {
    o << "[ i_idx = "; for(int i: b.i_idx) o << i << ", "; o << endl;
    o << "  j_idx = "; for(int j: b.j_idx) o << j << ", "; o << "]" << endl;
  }
  o << "END BLOCK LIST" << endl;
}

// ======================================================================
// For local low-rank approximation
// ======================================================================

/** Node supporting degenerate approximation by interpolation */
/* SAM_LISTING_BEGIN_Y */
template <int d>
class InterpNode: public CtNode<d> {
public:
  // Constructor from sequence of points; initializes \cob{$\VV$}
  InterpNode(const vector<Point<d>> _pts,size_t _q,int _dir=0):
    CtNode<d>(_pts,_dir),q(_q),sons{nullptr,nullptr},k(_pts.size()) 
  { initV(); } 
  virtual ~InterpNode(void) {} 
protected:
  // Initialization of matrix \cob{$\VV$}
  void initV(void);
public:
  const int q; // Rank, no of interpolation nodes 
  size_t k;    // Number of indices contained 
  MatrixXd V;  // low-rank factor \cob{$\VV\in\bbR^{k,q}$}
  InterpNode *sons[2]; // Pointers to sons (of type InterpNode!)
};
/* SAM_LISTING_END_Y */

/* SAM_LISTING_BEGIN_Z */
template <int d>
void InterpNode<d>::initV(void) { 
  static_assert(d == 1,"Implemented only for 1D");
  cerr << "Rank-" << q << " InterpNode: Initialization of V not implemented" << endl;
}

/* SAM_LISTING_END_Z */

template<int dim>
ostream &operator << (ostream &o,const InterpNode<dim> &node) {
  o << "## IPNode: " << node.pts.size()  << " points, " << node.getBBox() << ": ";
  for(const Point<dim> &pt: node.pts) {
    o << "["; for(int l=0;l<dim;l++) { o << pt.x[l] << ' '; } o << "]";
  }
  o << endl;
  if (node.sons[0]) o << *(node.sons[0]);
  if (node.sons[1]) o << *(node.sons[1]);
  return o;
}

/** Extended class for cluster trees for local low-rank approximation */
/* SAM_LISTING_BEGIN_E */
template <class Node>
class LowRankClusterTree: public ClusterTree<Node> {
public:
  // Idle constructor just setting rank argument q
  LowRankClusterTree(size_t _q):q(_q) {} 
  // Actual constructor taking a sequence of points
  void init(const vector<Point<Node::dim>> pts,size_t minpts = 1);
  virtual ~LowRankClusterTree(void) {}
protected:
  // factory method for relevant type of node takine rank argument
  virtual Node* createNode(const vector<Point<Node::dim>> pts,int dir) {
    return new Node(pts,q,dir); }
public:
  const size_t q; // rank of degenerate approximation on cluster boxes
};

template <class Node>
void LowRankClusterTree<Node>::init(const vector<Point<Node::dim>> pts,size_t minpts)
{ ClusterTree<Node>::init(pts,minpts); }
/* SAM_LISTING_END_E */

/** Type for far field cluster */ 
/* SAM_LISTING_BEGIN_F */
template <class Node,typename KERNEL>
class BiDirChebInterpBlock: public IndexBlock<Node> {
public:
  using kernel_t = KERNEL;
  BiDirChebInterpBlock(const Node &nx,const Node &ny,
		       kernel_t _G,size_t _q);
  virtual ~BiDirChebInterpBlock(void) {}
  // Invalid constructor throwing exception
  BiDirChebInterpBlock(const Node &nx,const Node &ny);
  
  const kernel_t G; // kernel function \cob{$\krn$}
  const int q; // No of interpolation nodes 
  MatrixXd C;  // \cob{$\VC\in\bbR^{q,q}$}   
};
/* SAM_LISTING_END_F */

template <class Node,typename KERNEL> BiDirChebInterpBlock<Node,KERNEL>::
BiDirChebInterpBlock(const Node &_nx,const Node &_ny,kernel_t _G,size_t _q):
  IndexBlock<Node>(_nx,_ny),G(_G),q(_q) 
{
  cerr << "BiDirChebInterpBlock: q = " << q << ": Initialization not yet implemented" << endl;
}

template <class Node,typename KERNEL> BiDirChebInterpBlock<Node,KERNEL>::
BiDirChebInterpBlock(const Node &_nx,const Node &_ny):
  IndexBlock<Node>(_nx,_ny),G(kernel_t()),q(0) 
{
  throw runtime_error("BiDirChebInterpBlock: no construction without kernel");
}

/** General type for generic near-field cluster pair */
/* SAM_LISTING_BEGIN_G */
template <class Node,typename KERNEL>
class NearFieldBlock: public IndexBlock<Node> {
public:
  using kernel_t = KERNEL;
  NearFieldBlock(const Node &nx,const Node &ny,kernel_t _G);
  virtual ~NearFieldBlock(void) {}
  // Invalid constructor throwing exception
  NearFieldBlock(const Node &nx,const Node &ny);
  
  const kernel_t G; // kernel function \cob{$\krn$}
  MatrixXd Mloc;   // local kernel collocation matrix
};
/* SAM_LISTING_END_G */

template <class Node,typename KERNEL> NearFieldBlock<Node,KERNEL>::
NearFieldBlock(const Node &_nx,const Node &_ny,kernel_t _G):
  IndexBlock<Node>(_nx,_ny),G(_G)
{
  cerr << "NearFieldBlock: Initialization not yet implemented" << endl;
}

template <class Node,typename KERNEL> NearFieldBlock<Node,KERNEL>::
NearFieldBlock(const Node &_nx,const Node &_ny):
  IndexBlock<Node>(_nx,_ny),G(kernel_t())
{
  throw runtime_error("Near Field Block: no construction without kernel");
}

/** Extended class for block partition, knowing low-rank compression */
/* SAM_LISTING_BEGIN_H */
template <class Node,typename FFB,typename NFB>
class HierMatBlockPartition: public BlockPartition<Node,FFB,NFB> {
public:
  using kernel_t = typename NFB::kernel_t;
  HierMatBlockPartition(const LowRankClusterTree<Node> &_xT,
			const LowRankClusterTree<Node> &_yT,
			kernel_t _G,size_t _q,
			double eta0 = 0.5):
    BlockPartition<Node,FFB,NFB>(_xT,_yT),G(_G),q(_q)
  { BlockPartition<Node,FFB,NFB>::init(eta0); }
  virtual ~HierMatBlockPartition(void) {}
protected:
  // Construct an instance of far-field block type
  virtual FFB makeFarFieldBlock(const Node &nx,const Node &ny) const
  {
    cout << "HierMatBlockPartition: makeFarFieldBlock" << endl;
    return FFB(nx,ny,G,q);
  }
  // Construct an instance of near-field block type
  virtual NFB makeNearFieldBlock(const Node &nx,const Node &ny) const
  {
    cout << "HierMatBlockPartition: makeNearFieldBlock" << endl;
    return NFB(nx,ny,G);
  }
public:
  const kernel_t G; // Reference to the kernel function
  const size_t q;  // degree+1 of interpolating polynomial
};
/* SAM_LISTING_END_H */

/** Class providing the kernel function */
struct Kernel {
  double operator () (double x,double y) {
    if (x != y) return -std::log(std::abs(x-y)); else return 0.0;
  }
};

/* SAM_LISTING_BEGIN_R */
template<int d>  // dimension \cob{$d$} as template argument
struct BasisFn {
  size_t idx;                   // Index of basis function
  Matrix<double,d,1> xmin,xmax; // Corners of bounding box
};
/* SAM_LISTING_END_R */

int main(int,char **) {
  const int d=1;
  size_t npts = 16;
  vector<Point<d>> pts;
  for(int n=0;n<npts;n++) {
    Point<d> p; p.idx = n; p.x[0] = (double)n;
    pts.push_back(p);
  }
  {
    ClusterTree<CtNode<d>> T; T.init(pts); 
    cout << T << endl;

    BlockPartition<CtNode<d>,IndexBlock<CtNode<d>>,IndexBlock<CtNode<d>> > bP(T,T);
    bP.init();
    cout << bP;
    cout << "Part I: Normal termination" << endl;
  }
  {
    size_t q = 7;
    Kernel G;
    LowRankClusterTree<InterpNode<d>> LRT(q); LRT.init(pts); 
    cout << "LowRandClusterTree:" << endl << LRT << endl;
    HierMatBlockPartition<InterpNode<d>,
			  BiDirChebInterpBlock<InterpNode<d>,Kernel>,
			  NearFieldBlock<InterpNode<d>,Kernel>> lrbP(LRT,LRT,G,q);
    cout << "Part II: Normal termination" << endl;
    
  }
  exit(0);
}

  
  
