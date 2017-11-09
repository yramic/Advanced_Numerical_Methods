/***********************************************************************
 *                                                                     *
 * Code for Course "Advanced Numerical Methods for CSE"                *
 * (Prof. Dr. R. Hiptmair)                                             * 
 * Author: R.H.                                                        *
 * Date: Nov 8, 2017                                                   *
 * (C) Seminar for Applied Mathematics, ETH Zurich                     *
 * This code can be freely used for non-commercial purposes as long    *
 * as this header is left intact.                                      *
 ***********************************************************************/

// General includes
#include <iostream>
#include <cmath>
#include <limits>
#include <vector>
#include <exception>
#include <Eigen/Dense>

using namespace std;
using namespace Eigen;

/** @brief Data structure for a collocation point
   A collocation point has an index and coordinates 
   @param t dimension of ambient space */
/* SAM_LISTING_BEGIN_1 */
template<int d>
struct Point {
  size_t idx;            // number of collocation point
  Matrix<double,d,1> x;  // coordinate vector
};
/* SAM_LISTING_END_1 */

/** @brief Data structure for bounding box
   a bounding box is defined by two position vectors */
/* SAM_LISTING_BEGIN_2 */
template<int d>
struct BBox {
  // Bounding box from sequence of points
  BBox(const vector<Point<d>> pts);
  // Corner points of bounding box
  Matrix<double,d,1> minc,maxc;
};
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
  // Constructors taking a sequence of points
  CtNode(const vector<Point<d>> _pts,int _dir=0):
    pts(_pts),sons{nullptr,nullptr},dir(_dir) {}
  template <typename PTIT> CtNode(PTIT first,PTIT last,int _dir=0):
    pts(first,last),sons{nullptr,nullptr},dir(_dir) {}
  // Destructor (also attempts to destroy sons!
  virtual ~CtNode(void) {
    if (sons[0]) delete sons[0];
    if (sons[1]) delete sons[1];
  }
  // Access to bounding box
  BBox<d> getBBox(void) const { return BBox<d>(pts); }
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
template <class Node,int d>
class ClusterTree {
public:
  // Constructor and destructor
  ClusterTree(const vector<Point<d>> pts,size_t minpts = 1);
  virtual ~ClusterTree(void) { delete root; }
  // Output of tree
  template <class _Node,int dim>
  friend ostream &operator << (ostream &o,const ClusterTree<_Node,dim> &T);
private:
  // Recursive construction
  void buildRec(Node *nptr,size_t minpts);
public:
  Node *root;
};
/* SAM_LISTING_END_5 */

/* SAM_LISTING_BEGIN_6 */
template<class Node,int d>
ClusterTree<Node,d>::ClusterTree(const vector<Point<d>> pts,size_t minpts) {
  if (!(root = new Node(pts))) throw(runtime_error("Cannot allocate root"));
  buildRec(root,minpts);
}
/* SAM_LISTING_END_6 */

template <class Node,int dim>
ostream &operator << (ostream &o,const ClusterTree<Node,dim> &T) {
  return o << "ROOT: " << *(T.root) << "END CLUSTER TREE" << endl;
}
  
  
/* SAM_LISTING_BEGIN_7 */
template<class Node,int d>
void ClusterTree<Node,d>::buildRec(Node *nptr,size_t minpts) {
  const size_t n = nptr->pts.size(); // Number of points
  if (n > minpts) {
    // Points have to be copied and sorted according to direction dir
    vector<Point<d>> tpts(nptr->pts);
    const int dir = (nptr->dir + 1)%d;
    sort(tpts.begin(),tpts.end(),
	 [dir] (const Point<d> &p1,const Point<d> &p2) -> bool {
	   return (bool)(p1.x[dir] < p2.x[dir]); });
    // Split point sequence and construct sons
    const size_t m = n/2; // integer arithmeric !
    auto first=tpts.cbegin();
    if (m > 0) {
      if (!(nptr->sons[0] = new Node(first,first+m,dir)))
	throw(runtime_error("Cannot allocate first son"));
      buildRec(nptr->sons[0],minpts); // recurse into first son
    } 
    if (!(nptr->sons[1] = new Node(first+m,tpts.cend(),dir)))
      throw(runtime_error("Cannot allocate second son"));
    buildRec(nptr->sons[1],minpts); // recurse into 2nd son
  }
}
/* SAM_LISTING_END_7 */

int main(int,char **) {
  const int d=1;
  size_t npts = 16;
  vector<Point<d>> pts;
  for(int n=0;n<npts;n++) {
    Point<d> p; p.idx = n; p.x[0] = (double)n;
    pts.push_back(p);
  }
  ClusterTree<CtNode<d>,d> T(pts);
  cout << T << endl;

  exit(0);
}

  
  
