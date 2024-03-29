/***********************************************************************
 *                                                                     *
 * Code for Course "Advanced Numerical Methods for CSE"                *
 * (Prof. Dr. R. Hiptmair)                                             *
 * Author: Ioannis Magkanaris                                          *
 * Date: 11/2017                                                       *
 * (C) Seminar for Applied Mathematics, ETH Zurich                     *
 * This code can be freely used for non-commercial purposes as long    *
 * as this header is left intact.                                      *
 ***********************************************************************/
#include "../../include/node.hpp"

#include <Eigen/SVD>
#include <iostream>

#include "../../include/cheby.hpp"
#include "../../include/point.hpp"

#define equal_clusters

// actual  constructor: creates the root of the Cluster Tree and the recursivly creates the leaves
Node::Node(std::vector<Point> Points, unsigned deg)
    : tl_child_(NULL),
      tr_child_(NULL),
      bl_child_(NULL),
      br_child_(NULL),
      deg_(deg),
      PPointsTree_(Points),
      CVc_node_(Eigen::VectorXd::Zero((deg + 1) * (deg + 1))) {
  setSons();
}

// destructor
Node::~Node() {
  if (tl_child_ != NULL) delete tl_child_;
  if (tr_child_ != NULL) delete tr_child_;
  if (bl_child_ != NULL) delete bl_child_;
  if (br_child_ != NULL) delete br_child_;
}

// calculate the rectangle defined by the points of the node
void Node::getRect() {
  /* SAM_LISTING_BEGIN_2 */
  double maxX, minX, maxY, minY;
  maxX = PPointsTree_.begin()->getX();
  minX = PPointsTree_.begin()->getX();
  maxY = PPointsTree_.begin()->getY();
  minY = PPointsTree_.begin()->getY();
  for (std::vector<Point>::iterator it = PPointsTree_.begin() + 1;
       it != PPointsTree_.end(); it++) {
    if (it->getX() > maxX) {
      maxX = it->getX();
    } else if (it->getX() < minX) {
      minX = it->getX();
    }
    if (it->getY() > maxY) {
      maxY = it->getY();
    } else if (it->getY() < minY) {
      minY = it->getY();
    }
  }
  x1_ = minX;
  x2_ = maxX;
  y1_ = minY;
  y2_ = maxY;
  Cheby cbx(x1_, x2_, deg_);
  Cheby cby(y1_, y2_, deg_);
  tkx_ = cbx.getNodes();  // Chebyshew nodes for x axis
  wkx_ = cbx.getWghts();  // weights of Lagrange polynomial for x axis
  tky_ = cby.getNodes();  // Chebyshew nodes for y axis
  wky_ = cby.getWghts();  // weights of Lagrange polynomial for y axis
                          /* SAM_LISTING_END_2 */
}

// build tree recursively
void Node::setSons() {
  if (!PPointsTree_.empty() &&
      PPointsTree_.size() >
          1) {  // if there are points in the PPointsTree vector of points then they are equaly divided into the node´s children
#ifdef equal_clusters
    /* SAM_LISTING_BEGIN_0 */
    auto checkX = [](Point a, Point b) -> bool { return a.getX() < b.getX(); };
    std::sort(PPointsTree_.begin(), PPointsTree_.end(), checkX);
    std::vector<Point>::iterator it;
    it =
        PPointsTree_.begin() +
        (PPointsTree_.size() + 1) /
            2;  // set iterator in the middle of the vector of points of this node
    std::vector<Point> l_points,
        r_points;  // now sort l\_points and r\_points based on their y-coordinates
    l_points.assign(PPointsTree_.begin(), it);
    r_points.assign(it, PPointsTree_.end());
    auto checkY = [](Point a, Point b) -> bool { return a.getY() < b.getY(); };
    std::sort(
        l_points.begin(), l_points.end(),
        checkY);  // sort left and right vectors into top and bottom based on the y-coordinates
    std::sort(r_points.begin(), r_points.end(), checkY);
    std::vector<Point> tl_PPoints, tr_PPoints, bl_PPoints,
        br_PPoints;  // creation of vectors of points of child nodes
    it = l_points.begin() + (l_points.size() + 1) / 2;
    tl_PPoints.assign(it, l_points.end());
    bl_PPoints.assign(l_points.begin(), it);
    it = r_points.begin() + (r_points.size() + 1) / 2;
    tr_PPoints.assign(it, r_points.end());
    br_PPoints.assign(r_points.begin(), it);

    if (!tl_PPoints.empty())
      tl_child_ = new Node(
          tl_PPoints,
          deg_);  // recursive construction of the Cluster Tree levels below root
    if (!tr_PPoints.empty()) tr_child_ = new Node(tr_PPoints, deg_);
    if (!bl_PPoints.empty()) bl_child_ = new Node(bl_PPoints, deg_);
    if (!br_PPoints.empty()) br_child_ = new Node(br_PPoints, deg_);
      /* SAM_LISTING_END_0 */
#endif
#ifdef inertia
    /* SAM_LISTING_BEGIN_1 */
    double sumX, sumY;
    for (int i = 0; i < PPointsTree_.size(); i++) {
      sumX += PPointsTree_[i].getX();
      sumY += PPointsTree_[i].getY();
    }
    double avgX = sumX / (double)PPointsTree_.size(),
           avgY = sumY / (double)PPointsTree_.size();
    Eigen::MatrixXd A(2, PPointsTree_.size());
    for (unsigned i = 0; i < PPointsTree_.size(); ++i) {
      A(0, i) = PPointsTree_[i].getX() - avgX;
      A(1, i) = PPointsTree_[i].getY() - avgY;
    }
    Eigen::MatrixXd M = A * A.transpose();
    Eigen::JacobiSVD<Eigen::MatrixXd> svdOfM(
        M,
        Eigen::
            ComputeThinV);  // 'M' is square, so ComputeFullV and ComputeThinV are the same
    Eigen::MatrixXd V = svdOfM.matrixV();
    auto y = [](double x, double x1, double y1, double avgX,
                double avgY) -> double { return avgY + (x - avgX) * y1 / x1; };
    std::vector<Point> top_PPoints, bottom_PPoints, left_PPoints,
        right_PPoints;  // creation of vectors of points of child nodes
    for (unsigned i = 0; i < PPointsTree_.size(); ++i) {
      double y1 = y(
          PPointsTree_[i].getX(), V(0, 0), V(1, 0), avgX,
          avgY);  // y value of the line that is defined by the point {avgX,avgY} and the vector corresponding to the biggest eigen value
      double y2 = y(
          PPointsTree_[i].getX(), V(0, 1), V(1, 1), avgX,
          avgY);  // y value of the line that is defined by the point {avgX,avgY} and the vector corresponding to the second biggest eigen value
      if (y2 <= PPointsTree_[i].getY()) {
        if (y1 <= PPointsTree_[i].getY()) {
          top_PPoints.push_back(PPointsTree_[i]);
        } else {
          right_PPoints.push_back(PPointsTree_[i]);
        }
      } else {
        if (y1 <= PPointsTree_[i].getY()) {
          left_PPoints.push_back(PPointsTree_[i]);
        } else {
          bottom_PPoints.push_back(PPointsTree_[i]);
        }
      }
    }

    if (!top_PPoints.empty())
      tl_child_ = new Node(
          top_PPoints,
          deg_);  // recursive construction of the Cluster Tree levels below root
    if (!right_PPoints.empty()) tr_child_ = new Node(right_PPoints, deg_);
    if (!bottom_PPoints.empty()) bl_child_ = new Node(bottom_PPoints, deg_);
    if (!left_PPoints.empty()) br_child_ = new Node(left_PPoints, deg_);
      /* SAM_LISTING_END_1 */
#endif
    getRect();  // calculate the rectangle defined by the points of the node
    // fix for the rectangle if it is a segment
    if (std::abs(x1_ - x2_) < 10 * std::numeric_limits<double>::epsilon()) {
      x2_++;
    }
    if (std::abs(y1_ - y2_) < 10 * std::numeric_limits<double>::epsilon()) {
      y2_++;
    }
    Cheby cbx(x1_, x2_,
              deg_);  // Chebvchev interpolation on the edges of the rectangle
    Cheby cby(y1_, y2_, deg_);
    tkx_ = cbx.getNodes();  // Chebyshew nodes for x axis
    wkx_ = cbx.getWghts();  // weights of Lagrange polynomial for x axis
    tky_ = cby.getNodes();  // Chebyshew nodes for y axis
    wky_ = cby.getWghts();  // weights of Lagrange polynomial for y axis
  }
}

// compute V-matrix of node
unsigned Node::setV() {
  /* SAM_LISTING_BEGIN_3 */
  unsigned ppts = PPointsTree_.size();
  auto checkID = [](Point a, Point b) -> bool { return a.getId() < b.getId(); };
  std::sort(PPointsTree_.begin(), PPointsTree_.end(), checkID);
  V_node_ = Eigen::MatrixXd::Constant(ppts, (deg_ + 1) * (deg_ + 1), 1);
  for (unsigned i = 0; i <= ppts - 1;
       ++i) {  // calculation of Vx combined with Vy
    for (unsigned j1 = 0; j1 <= deg_; ++j1) {
      for (unsigned k1 = 0; k1 < j1; ++k1) {
        for (unsigned j2 = 0; j2 <= deg_; ++j2) {
          V_node_(i, j1 * (deg_ + 1) + j2) *=
              (PPointsTree_[i].getX() - tkx_[k1]);
        }
      }
      // Skip "k1 == j1"
      for (unsigned k1 = j1 + 1; k1 <= deg_; ++k1) {
        for (unsigned j2 = 0; j2 <= deg_; ++j2) {
          V_node_(i, j1 * (deg_ + 1) + j2) *=
              (PPointsTree_[i].getX() - tkx_[k1]);
        }
      }
      for (unsigned j2 = 0; j2 <= deg_; ++j2) {
        for (unsigned k2 = 0; k2 < j2; ++k2) {
          V_node_(i, j1 * (deg_ + 1) + j2) *=
              (PPointsTree_[i].getY() - tky_[k2]);
        }
        // Skip "k2 == j2"
        for (unsigned k2 = j2 + 1; k2 <= deg_; ++k2) {
          V_node_(i, j1 * (deg_ + 1) + j2) *=
              (PPointsTree_[i].getY() - tky_[k2]);
        }
        V_node_(i, j1 * (deg_ + 1) + j2) *= wkx_[j1] * wky_[j2];
      }
    }
  }

  // Alternate way of computing V matrix
  /*
    Eigen::MatrixXd VnodeX = Eigen::MatrixXd::Ones(ppts, (deg_+1));
    Eigen::MatrixXd VnodeY = Eigen::MatrixXd::Ones(ppts, (deg_+1));
    for(unsigned i=0; i<=ppts-1; ++i) {
        for(unsigned j=0; j<=deg_; ++j) {
            for(unsigned k=0; k<j; ++k) {
                VnodeX(i,j) *= PPointsTree_[i].getX() - tkx_[k];
            }
            // Skip "k == j"
            for(unsigned k=j+1; k<=deg_; ++k) {
                VnodeX(i,j) *= PPointsTree_[i].getX() - tkx_[k];
            }
            VnodeX(i,j) *= wkx_(j);
        }
    }
    for(unsigned i=0; i<=ppts-1; ++i) {
        for(unsigned j=0; j<=deg_; ++j) {
            for(unsigned k=0; k<j; ++k) {
                VnodeY(i,j) *= PPointsTree_[i].getY() - tky_[k];
            }
            // Skip "k == j"
            for(unsigned k=j+1; k<=deg_; ++k) {
                VnodeY(i,j) *= PPointsTree_[i].getY() - tky_[k];
            }
            VnodeY(i,j) *= wky_(j);
        }
    }

    Eigen::MatrixXd V_node_new(ppts, (deg_+1)*(deg_+1));
    for(unsigned i=0; i<=ppts-1; ++i) {
        for(unsigned j=0; j<=deg_; ++j) {
            V_node_new.block(i, j*(deg_+1), 1, deg_+1) = VnodeX(i,j) * VnodeY.row(i);
        }
    }
    V_node_ = V_node_new;
*/
  /* SAM_LISTING_END_3 */
  return ppts * (deg_ + 2) * (deg_ + 1) /
         2;  // return no. of 'operations' performed
}

// compute V*c restricted to node indices
unsigned Node::setVc(const Eigen::VectorXd& c) {
  int n = PPointsTree_.size();
  Eigen::VectorXd c_seg = Eigen::VectorXd::Zero(n);
  for (int i = 0; i < n; i++) {  // get only the part of vector c needed
    c_seg[i] = c(PPointsTree_[i].getId());
  }
  Vc_node_ = V_node_.transpose() * c_seg;  // Vc matrix calculation
  return V_node_.rows() *
         V_node_.cols();  // return no. of 'operations' performed
}

// update C*V*c restricted to node indices
unsigned Node::setCVc(const Eigen::VectorXd& CVc) {
  CVc_node_ += CVc;
  return CVc_node_.size();  // return no. of 'operations' performed
}

void Node::printree(int n) {
  std::cout << "Node " << n << std::endl;
  for (std::vector<Point>::iterator it = PPointsTree_.begin();
       it != PPointsTree_.end(); it++) {
    std::cout << it->getId() << " ";
  }
  std::cout << std::endl;
  if (tl_child_ != NULL) {
    std::cout << "tl " << n << " ";
    tl_child_->printree(n + 1);
  }
  if (tr_child_ != NULL) {
    std::cout << "tr " << n << " ";
    tr_child_->printree(n + 1);
  }
  if (bl_child_ != NULL) {
    std::cout << "bl " << n << " ";
    bl_child_->printree(n + 1);
  }
  if (br_child_ != NULL) {
    std::cout << "br " << n << " ";
    br_child_->printree(n + 1);
  }
}
