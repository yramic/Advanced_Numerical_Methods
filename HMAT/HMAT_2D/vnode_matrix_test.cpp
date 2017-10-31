#include "include/point.hpp"
#include "include/ctree.hpp"
#include "include/node.hpp"
#include "include/cheby.hpp"
#include "include/low_rank_app.hpp"
#include <iostream>

int main() {
    unsigned n=16;  // number of points
    std::vector<int> x = {24,22,73,63,14,17,39,99,83,41,40,4,83,30,65,23}, y = {24,55,87,91,30,1,9,28,67,85,10,84,57,72,86,56}, v = {66,63,27,46,83,58,46,8,95,57,2,79,34,21,64,95};
    std::vector<Point> PPoints; // initalizing Polygon Points properties
    PPoints.reserve(n);
    for (int i=0; i<n; i++){
        Point p;
        p.setId(i);
        p.setX(x[i]);
        p.setY(y[i]);
        p.setV(v[i]);
        PPoints.push_back(p);
    }
    unsigned deg = 2;
    Node t(PPoints,deg);
    double eta;
    PolynomialKernel P;
    LowRankApp HMat(&P, PPoints, n, eta, deg);
    Eigen::MatrixXd V = HMat.setV_node(&t,deg);
    // alternate calculation of the V matrix of this Node
    double x1,x2,y1,y2;                 // construction of Bounding Box of this Node
    x1 = PPoints.begin()->getX();
    x2 = PPoints.begin()->getX();
    y1 = PPoints.begin()->getY();
    y2 = PPoints.begin()->getY();
    for (std::vector<Point>::iterator it=PPoints.begin(); it!=PPoints.end(); it++) {
        if (it->getX()>x2) {
            x2 = it->getX();
        }
        else if (it->getX()<x1) {
            x1 = it->getX();
        }
        if (it->getY()>y2) {
            y2 = it->getY();
        }
        else if (it->getY()<y1) {
            y1 = it->getY();
        }
    }
    double x1_b_ = x1;
    double x2_b_ = x2;
    double y1_b_ = y1;
    double y2_b_ = y2;
    // fix for points of a bbox being a segment
    if(std::abs(x1-x2)<10*std::numeric_limits<double>::epsilon()){
        x2_b_++;
    }
    if(std::abs(y1-y2)<10*std::numeric_limits<double>::epsilon()){
        y2_b_++;
    }
    // calculate Vnode for x axis and then for y axis
    Cheby cbx(x1_b_, x2_b_, deg);
    Cheby cby(y1_b_, y2_b_, deg);
    Eigen::VectorXd tkx = cbx.getNodes(); // Chebyshew nodes for x axis
    Eigen::VectorXd wkx = cbx.getWghts(); // weights of Lagrange polynomial for x axis
    Eigen::VectorXd tky = cby.getNodes(); // Chebyshew nodes for y axis
    Eigen::VectorXd wky = cby.getWghts(); // weights of Lagrange polynomial for y axis
    Eigen::MatrixXd VnodeX = Eigen::MatrixXd::Constant(n, (deg+1), 1);
    Eigen::MatrixXd VnodeY = Eigen::MatrixXd::Constant(n, (deg+1), 1);
    for(unsigned i=0; i<=n-1; ++i) {
        for(unsigned j=0; j<=deg; ++j) {
            for(unsigned k=0; k<j; ++k) {
                VnodeX(i,j) *= PPoints[i].getX() - tkx[k];
            }
            // Skip "k == j"
            for(unsigned k=j+1; k<=deg; ++k) {
                VnodeX(i,j) *= PPoints[i].getX() - tkx[k];
            }
            VnodeX(i,j) *= wkx(j);
        }
    }
    for(unsigned i=0; i<=n-1; ++i) {
        for(unsigned j=0; j<=deg; ++j) {
            for(unsigned k=0; k<j; ++k) {
                VnodeY(i,j) *= PPoints[i].getY() - tky[k];
            }
            // Skip "k == j"
            for(unsigned k=j+1; k<=deg; ++k) {
                VnodeY(i,j) *= PPoints[i].getY() - tky[k];
            }
            VnodeY(i,j) *= wky(j);
        }
    }
    // calculate the V node
    Eigen::MatrixXd V_node_new(n, (deg+1)*(deg+1));
    for(unsigned i=0; i<=n-1; ++i) {
        for(unsigned j=0; j<=deg; ++j) {
            V_node_new.block(i, j*(deg+1), 1, deg+1) = VnodeX(i,j) * VnodeY.row(i);
        }
    }
    // check if the two ways of calculation produce the same V matrix
    for(int i = 0; i<n; i++){
        for(int j = 0; j<deg*deg; j++){
            if(std::abs(V_node_new(i,j) - V(i,j))>10*std::numeric_limits<double>::epsilon()) {
                std::cout << "Wrong" << std::endl;
                return 0;
            }
        }
    }
    std::cout << "Correct" << std::endl;
    return 0;
}
