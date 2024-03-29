/***********************************************************************
 *                                                                     *
 * Code for Course "Advanced Numerical Methods for CSE"                *
 * (Prof. Dr. R. Hiptmair)                                             *
 * Author: Daniele Casati                                              *
 * Date: 11/2017                                                       *
 * (C) Seminar for Applied Mathematics, ETH Zurich                     *
 * This code can be freely used for non-commercial purposes as long    *
 * as this header is left intact.                                      *
 ***********************************************************************/
#include "../include/low_rank_app.hpp"

#include <Eigen/Dense>
#include <chrono>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>

#include "../include/block_cluster.hpp"
#include "../include/block_nearf.hpp"
#include "../include/ctree.hpp"
#include "../include/kernel.hpp"
#include "../include/node.hpp"
#include "../include/point.hpp"

// constructor
template <>
LowRankApp<BlockCluster, Node>::LowRankApp(Kernel *kernel,
                                           const std::vector<Point> &GPoints,
                                           double eta, unsigned deg)
    : kernel_(kernel),
      GPoints_(GPoints),
      HP_(GPoints, eta, deg),
      deg_(deg),
      nops_(0),
      debug_(false) {
  // Recursive hierarchical partitioning to get Near/Far field boxes
  HP_.setNearFar();
}

// constructor
template <>
LowRankApp<BlockCluster, Node>::LowRankApp(Kernel *kernel,
                                           const std::vector<Point> &GPoints,
                                           double eta, unsigned deg,
                                           const std::string &filename)
    : kernel_(kernel),
      GPoints_(GPoints),
      HP_(GPoints, eta, deg),
      deg_(deg),
      nops_(0),
      debug_(true),
      myfile_(filename.c_str(), std::ios::app) {
  // Recursive hierarchical partitioning to get Near/Far field boxes
  HP_.setNearFar();
}

// pre-processing: initialize matrix V and vector Vc for all far field nodes
// do these steps once for each node, not every time the node appears in a pair
template <>
void LowRankApp<BlockCluster, Node>::preProcess(std::vector<Node *> ff_v_x,
                                                std::vector<Node *> ff_v_y,
                                                const Eigen::VectorXd &c) {
  for (auto &xnode : ff_v_x) {  // iterate for all the far field xnodes
    nops_ += xnode->setV();
  }
  for (auto &ynode : ff_v_y) {  // iterate for all the far field ynodes
    nops_ += ynode->setV();
    nops_ += ynode->setVc(c);
  }
}

// block-processing: compute vector CVc for all far field pairs and store it
// into xnode all vectors CVc of an xnode can already be summed together
template <>
void LowRankApp<BlockCluster, Node>::blockProcess(
    std::vector<BlockCluster *> ff_v) {
  for (auto &pair : ff_v) {  // iterate for all the pairs of far field nodes
    nops_ += pair->setMatrix(
        kernel_);  // here because needed for each pair of nodes,
                   // cannot be moved to pre-processing
    nops_ += pair->setCVc();
  }
}

// Sets CVc for the pair to 0
template <>
void LowRankApp<BlockCluster, Node>::blockProcessClear(
    std::vector<BlockCluster *> ff_v) {
  for (auto &pair : ff_v) {  // iterate for all the pairs of far field nodes
    pair->getXNode()->resetCVc();
  }
}

// debug-processing: compute approximate matrix VCV for all far field pairs,
// orresponding exact block C, and the error between them
template <>
void LowRankApp<BlockCluster, Node>::debugProcess(
    std::vector<BlockCluster *> ff_v) {
  double error_Frobenius = 0., error_max = 0.;
  for (auto &pair : ff_v) {  // iterate for all the pairs of far field nodes
    Eigen::MatrixXd block_approx = pair->getVCV();
    BlockNearF tmp(pair->getXNode(), pair->getYNode());
    tmp.setMatrix(kernel_);
    Eigen::MatrixXd block_exact = tmp.getMatrix();
    Eigen::MatrixXd diff_blocks = block_exact - block_approx;
    error_Frobenius +=
        diff_blocks.norm() / std::sqrt(block_exact.rows() * block_exact.cols());
    double error_tmp = diff_blocks.cwiseAbs().maxCoeff();
    if (error_max < error_tmp) {
      error_max = error_tmp;
    }
  }

  myfile_ << "error_Frobenius, " << GPoints_.size() << ", "
          << std::setprecision(10) << error_Frobenius / ff_v.size()
          << std::endl;  // average error w.r.t. all blocks
  myfile_ << "error_max, " << GPoints_.size() << ", " << std::setprecision(10)
          << error_max << std::endl;
}

// post-processing: compute vector Vx*CVc for all far field xnodes and add it to
// vector f in the right place
template <>
void LowRankApp<BlockCluster, Node>::postProcess(std::vector<Node *> ff_v_x,
                                                 Eigen::VectorXd &f) {
  for (auto &xnode : ff_v_x) {  // iterate for all the far field xnodes
    Eigen::VectorXd CVc = xnode->getCVc_Node();
    Eigen::MatrixXd Vx = xnode->getV_Node();
    Eigen::VectorXd f_seg = Vx * CVc;
    nops_ += Vx.rows() * Vx.cols();
    for (int i = 0; i < xnode->getPoints().size(); i++) {
      f[xnode->getPoints()[i].getId()] +=
          f_seg[i];  // add contribution of far field to ``f''
    }
  }
}

// count far-field ynodes contributing to each row of the approximate low-rank
// matrix
template <>
void LowRankApp<BlockCluster, Node>::calc_numb_approx_per_row(
    std::vector<BlockCluster *> ff_v, Eigen::VectorXd &f_approx_ff_contr) {
  for (auto &pair : ff_v) {  // iterate for all the pairs of far field nodes
    Node *xnode = pair->getXNode();
    Node *ynode = pair->getYNode();
    for (int i = 0; i < xnode->getPoints().size(); i++) {
      f_approx_ff_contr(xnode->getPoints()[i].getId()) +=
          ynode->getPoints().size();
    }
  }
}

// compute far-field contribution
template <>
void LowRankApp<BlockCluster, Node>::ff_contribution(
    std::vector<BlockCluster *> ff_v, std::vector<Node *> ff_v_x,
    std::vector<Node *> ff_v_y, const Eigen::VectorXd &c, Eigen::VectorXd &f,
    Eigen::VectorXd &f_approx_ff_contr) {
  preProcess(ff_v_x, ff_v_y, c);
  blockProcess(ff_v);
  if (debug_) {
    debugProcess(ff_v);
  }
  postProcess(ff_v_x, f);
  calc_numb_approx_per_row(ff_v, f_approx_ff_contr);
}

// compute near-field contribution
template <>
void LowRankApp<BlockCluster, Node>::nf_contribution(
    std::vector<BlockNearF *> nf_v, const Eigen::VectorXd &c,
    Eigen::VectorXd &f, Eigen::VectorXd &f_approx_nf_contr) {
  for (auto &pair : nf_v) {  // iterate for all the near field xnodes
    Node *xnode = pair->getXNode();
    Node *ynode = pair->getYNode();
    nops_ += pair->setMatrix(kernel_);
    Eigen::MatrixXd C = pair->getMatrix();
    for (int i = 0; i < xnode->getPoints().size(); i++) {
      for (int j = 0; j < ynode->getPoints().size(); j++) {
        f(xnode->getPoints()[i].getId()) +=
            C(i, j) * c(ynode->getPoints()[j]
                            .getId());  // add near field contribution to ``f''
        // The contributions of all near-field pairs involving 'xnode' can first
        // be summed ('blockProcess') and only then positioned in the right
        // entry of 'f' ('postProcess'), similarly to the near-field vector.
        ++f_approx_nf_contr(xnode->getPoints()[i].getId());
        nops_ += C.rows() * C.cols();
      }
    }
  }
}

// approximate matrix-vector multiplication
template <>
Eigen::VectorXd LowRankApp<BlockCluster, Node>::mvProd(
    const Eigen::VectorXd &c) {
  blockProcessClear(HP_.getFF());
  // Setting V and Vc for far field nodes
  preProcess(HP_.getFFxnds(), HP_.getFFynds(), c);
  // Setting CVc for far field blocks
  blockProcess(HP_.getFF());

  Eigen::VectorXd f_approx_ff_contr = Eigen::VectorXd::Zero(c.size());
  Eigen::VectorXd f_approx_nf_contr = Eigen::VectorXd::Zero(c.size());

  // compute far field contribution
  Eigen::VectorXd f_approx = Eigen::VectorXd::Zero(c.size());
  //auto start = std::chrono::high_resolution_clock::now();
  // Actual multiplication
  postProcess(HP_.getFFxnds(), f_approx);
  // std::cout << "Far Field Contribution for each row" << std::endl;
  // std::cout << f_approx_ff_contr << std::endl;

  // compute near-field contribution
  nf_contribution(HP_.getNF(), c, f_approx, f_approx_nf_contr);
  // std::cout << "Near Field Contribution for each row" << std::endl;
  // std::cout << f_approx_nf_contr << std::endl;
  //auto end = std::chrono::high_resolution_clock::now();
  //auto time_diff =
  //    std::chrono::duration_cast<std::chrono::microseconds>(end - start);
  //f_approx(0) = time_diff.count();
  // std::cout << "Near Field Nodes: " << near << " Far Field Nodes: " << far <<
  // std::endl; std::cout << "Near Field Nodes: " << (double)near/(near+far)*100.
  // << "% " << "Far Field Nodes: " << (double)far/(near+far)*100. << "%" <<
  // std::endl;

  // std::cout << "Number of matrix operations performed for low-rank
  // approximation: " << nops_ << std::endl;

  return f_approx;
}
