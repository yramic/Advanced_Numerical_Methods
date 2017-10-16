#include "../include/node.hpp"
#include "../include/cheby.hpp"
#include <Eigen/Dense>
#include <iostream>
#include "../include/point.hpp"
#include <limits>

// actual  constructor: creates the root of the Cluster Tree and the recursivly creates the leaves
Node::Node(std::vector<Point> PPointsTree):
    l_child_(NULL), r_child_(NULL), l_ind_(0), r_ind_(0), tl_child_(NULL), tr_child_(NULL), bl_child_(NULL), br_child_(NULL), PPointsTree_(PPointsTree), near_f_(vector_t()), far_f_(vector_t())//, x1_(-(std::numeric_limits<double>::max())), x2_(std::numeric_limits<double>::max()), y1_(-(std::numeric_limits<double>::max())), y2_(std::numeric_limits<double>::max())
{
    setLeaves();
}

// recursive constructor for the leaves of the Cluster Tree
Node::Node(std::vector<Point> PPointsTree, double x1, double x2, double y1, double y2):
    l_child_(NULL), r_child_(NULL), l_ind_(0), r_ind_(0), tl_child_(NULL), tr_child_(NULL), bl_child_(NULL), br_child_(NULL), PPointsTree_(PPointsTree), near_f_(vector_t()), far_f_(vector_t()), x1_(x1), x2_(x2), y1_(y1), y2_(y2)
{
    setLeaves(x1,x2,y1,y2);
}

// actual constructor: adds a tree below the node if left_index != right_index
Node::Node(unsigned l_ind, unsigned r_ind):
    tl_child_(NULL), tr_child_(NULL), bl_child_(NULL), br_child_(NULL), PPointsTree_(std::vector<Point>()), l_child_(NULL), r_child_(NULL), l_ind_(l_ind), r_ind_(r_ind), near_f_(vector_t()), far_f_(vector_t())
{
    setLeaves();
}


// destructor
Node::~Node()
{
    if((l_child_ == NULL) && (r_child_ == NULL))
        std::cout << "leaves destroyed" << std::endl;
    if(l_child_ != NULL)
        delete l_child_;
    if(r_child_ != NULL)
        delete r_child_;
}


// build tree recursively
void Node::setLeaves()
{
    if(r_ind_ - l_ind_ > 0) {
        l_child_ = new Node(l_ind_, (l_ind_+r_ind_)/2);
        r_child_ = new Node((l_ind_+r_ind_)/2 + 1, r_ind_);
    }
    if(!PPointsTree_.empty() && PPointsTree_.size()>1){ // in the beggining the space must be divided based on the max coordinates of the points
        double medX,medY,maxX,minX,maxY,minY;
        maxX = PPointsTree_.begin()->getX();
        minX = PPointsTree_.begin()->getX();
        maxY = PPointsTree_.begin()->getY();
        minY = PPointsTree_.begin()->getY();
        for (std::vector<Point>::iterator it=PPointsTree_.begin()+1; it!=PPointsTree_.end(); it++){
            if (it->getX()>maxX) {
                maxX = it->getX();
            }
            else if (it->getX()<minX) {
                minX = it->getX();
            }
            if (it->getY()>maxY) {
                maxY = it->getY();
            }
            else if (it->getY()<minY) {
                minY = it->getY();
            }
        }
        x1_ = x1_b_ = minX;
        x2_ = x2_b_ = maxX;
        y1_ = y1_b_ = minY;
        y2_ = y2_b_ = maxY;

        //std::cout << sumX << " " << sumY << std::endl << std::flush;
        medX = minX+(maxX-minX)/2;                      // the space is devided in 4 quadrants and the point of each quadrant is saved in a vector on the coresponding Node
        medY = minY+(maxY-minY)/2;
        std::vector<Point> tl_PPoints, tr_PPoints, bl_PPoints, br_PPoints;
        for (std::vector<Point>::iterator it=PPointsTree_.begin(); it!=PPointsTree_.end(); it++){
            if (it->getX()<=medX && it->getY()<=medY) {
                bl_PPoints.push_back(*it);
            }
            else if (it->getX()>medX && it->getY()<medY) {
                br_PPoints.push_back(*it);
            }
            else if (it->getX()<=medX && it->getY()>=medY) {
                tl_PPoints.push_back(*it);
            }
            else if (it->getX()>medX && it->getY()>medY) {
                tr_PPoints.push_back(*it);
            }
            else {
                bl_PPoints.push_back(*it);
            }
            //std::cout << it->getId() << ' ';
        }
            //std::cout << x1_ << " " << x2_ << " " << y1_ << " " << y2_ << std::endl << std::flush;
            //std::cout << std::endl << std::flush;
            if (!tl_PPoints.empty()) tl_child_ = new Node(tl_PPoints,minX,medX,medY,maxY);   // recursive construction of the Cluster Tree levels below root
            if (!tr_PPoints.empty()) tr_child_ = new Node(tr_PPoints,medX,maxX,medY,maxY);
            if (!bl_PPoints.empty()) bl_child_ = new Node(bl_PPoints,minX,medX,minY,medY);
            if (!br_PPoints.empty()) br_child_ = new Node(br_PPoints,medX,maxX,minY,medY);
    }
}

void Node::setLeaves(double x1, double x2, double y1, double y2)
{
    if(r_ind_ - l_ind_ > 0) {
        l_child_ = new Node(l_ind_, (l_ind_+r_ind_)/2);
        r_child_ = new Node((l_ind_+r_ind_)/2 + 1, r_ind_);
        //std::cout << "wtf" << std::endl << std::flush;
    }
    if(!PPointsTree_.empty() && PPointsTree_.size()>1){ // below the root level the space [x1,x2]x[y1,y2] must be divided in quadrants recursivly
        //std::cout << sumX << " " << sumY << std::endl << std::flush;
        int medX = x1+(x2-x1)/2;
        int medY = y1+(y2-y1)/2;
        std::vector<Point> tl_PPoints, tr_PPoints, bl_PPoints, br_PPoints;
        for (std::vector<Point>::iterator it=PPointsTree_.begin(); it!=PPointsTree_.end(); it++){
            if (it->getX()<=medX && it->getY()<=medY) {
                bl_PPoints.push_back(*it);
            }
            else if (it->getX()>medX && it->getY()<medY) {
                br_PPoints.push_back(*it);
            }
            else if (it->getX()<=medX && it->getY()>=medY) {
                tl_PPoints.push_back(*it);
            }
            else if (it->getX()>medX && it->getY()>medY) {
                tr_PPoints.push_back(*it);
            }
            else {
                bl_PPoints.push_back(*it);
            }
            //std::cout << it->getId() << ' ';
        }
            //std::cout << x1_ << " " << x2_ << " " << y1_ << " " << y2_ << std::endl << std::flush;
            //std::cout << std::endl << std::flush;
            if (!tl_PPoints.empty()) tl_child_ = new Node(tl_PPoints,x1,medX,medY,y2);
            if (!tr_PPoints.empty()) tr_child_ = new Node(tr_PPoints,medX,x2,medY,y2);
            if (!bl_PPoints.empty()) bl_child_ = new Node(bl_PPoints,x1,medX,y1,medY);
            if (!br_PPoints.empty()) br_child_ = new Node(br_PPoints,medX,x2,y1,medY);
    }
}


// compute V-matrix of node
void Node::setV_node(const std::vector<Point> &t, unsigned deg) //tt==PPointsTree??
{
    //std::cout << "setV_Node test" << std::endl;
    /*if(r_ind_ - l_ind_ > 0) {

        double xmin = x[l_ind_];
        double xmax = x[r_ind_];
        Cheby cb(xmin, xmax, deg);
        Eigen::VectorXd tk = cb.getNodes(); // Chebyshew nodes
        Eigen::VectorXd wk = cb.getWghts(); // weights of Lagrange polynomial
        V_node_ = Eigen::MatrixXd::Constant(r_ind_-l_ind_+1, deg+1, 1);

        for(unsigned i=0; i<=r_ind_-l_ind_; ++i) {
            for(unsigned j=0; j<=deg; ++j) {
                for(unsigned k=0; k<j; ++k) {
                    V_node_(i,j) *= x[i+l_ind_] - tk[k];
                }
                // Skip "k == j"
                for(unsigned k=j+1; k<=deg; ++k) {
                    V_node_(i,j) *= x[i+l_ind_] - tk[k];
                }
                V_node_(i,j) *= wk(j);
            }
        }
    }*/
    if (PPointsTree_.size() > 1) {
        double x1,x2,y1,y2;                 // construction of Bounding Box of this Node
        x1 = PPointsTree_.begin()->getX();
        x2 = PPointsTree_.begin()->getX();
        y1 = PPointsTree_.begin()->getY();
        y2 = PPointsTree_.begin()->getY();
        for (std::vector<Point>::iterator it=PPointsTree_.begin(); it!=PPointsTree_.end(); it++) {
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
        x1_b_ = x1;
        x2_b_ = x2;
        y1_b_ = y1;
        y2_b_ = y2;
        Cheby cbx(x1, x2, deg);
        Cheby cby(y1, y2, deg);
        Eigen::VectorXd tkx = cbx.getNodes(); // Chebyshew nodes for x axis
        //std:: cout << "tkx" << tkx << std::endl;
        Eigen::VectorXd wkx = cbx.getWghts(); // weights of Lagrange polynomial for x axis
        //std:: cout << "wkx" << wkx << std::endl;
        Eigen::VectorXd tky = cby.getNodes(); // Chebyshew nodes for y axis
        //std:: cout << "tky" << tky << std::endl;
        Eigen::VectorXd wky = cby.getWghts(); // weights of Lagrange polynomial for y axis
        //std:: cout << "wky" << wky << std::endl;
        int ppts = PPointsTree_.size();
        V_node_ = Eigen::MatrixXd::Constant(ppts, (deg+1)*(deg+1), 1);

        for(unsigned i=0; i<=ppts-1; ++i) { // calculation of Vx combined with Vy
            for(unsigned j1=0; j1<=deg; ++j1) {
                for(unsigned k1=0; k1<j1; ++k1) {
                    for(unsigned j2=0; j2<=deg; ++j2) {
                        V_node_(i,j1*(deg+1) + j2) *= (PPointsTree_[i].getX() - tkx[k1]);
                        for(unsigned k2=0; k2<j2; ++k2) {
                            V_node_(i,j1*(deg+1) + j2) *= (PPointsTree_[i].getY() - tky[k2]);
                        }
                        // Skip "k2 == j2"
                        for(unsigned k2=j2+1; k2<=deg; ++k2) {
                            V_node_(i,j1*(deg+1) + j2) *= (PPointsTree_[i].getY() - tky[k2]);
                        }
                        V_node_(i,j1*(deg+1) + j2) *= wkx(j1)*wky(j2);
                    }
                }
                // Skip "k1 == j1"
                for(unsigned k1=j1+1; k1<=deg; ++k1) {
                    for(unsigned j2=0; j2<=deg; ++j2) {
                        V_node_(i,j1*(deg+1) + j2) *= (PPointsTree_[i].getX() - tkx[k1]);
                        for(unsigned k2=0; k2<j2; ++k2) {
                            V_node_(i,j1*(deg+1) + j2) *= (PPointsTree_[i].getY() - tky[k2]);
                        }
                        // Skip "k2 == j2"
                        for(unsigned k2=j2+1; k2<=deg; ++k2) {
                            V_node_(i,j1*(deg+1) + j2) *= (PPointsTree_[i].getY() - tky[k2]);
                        }
                        V_node_(i,j1*(deg+1) + j2) *= wkx(j1)*wky(j2);
                    }
                }
            }
        }
        /*Eigen::MatrixXd VnodeX = Eigen::MatrixXd::Constant(ppts, (deg+1), 1);
        Eigen::MatrixXd VnodeY = Eigen::MatrixXd::Constant(ppts, (deg+1), 1);
        for(unsigned i=0; i<=ppts-1; ++i) {
            for(unsigned j=0; j<=deg; ++j) {
                for(unsigned k=0; k<j; ++k) {
                    VnodeX(i,j) *= PPointsTree_[i].getX() - tkx[k];
                }
                // Skip "k == j"
                for(unsigned k=j+1; k<=deg; ++k) {
                    VnodeX(i,j) *= PPointsTree_[i].getX() - tkx[k];
                }
                VnodeX(i,j) *= wkx(j);
            }
        }
        for(unsigned i=0; i<=ppts-1; ++i) {
            for(unsigned j=0; j<=deg; ++j) {
                for(unsigned k=0; k<j; ++k) {
                    VnodeY(i,j) *= PPointsTree_[i].getY() - tky[k];
                }
                // Skip "k == j"
                for(unsigned k=j+1; k<=deg; ++k) {
                    VnodeY(i,j) *= PPointsTree_[i].getY() - tky[k];
                }
                VnodeY(i,j) *= wky(j);
            }
        }
        std::cout << "VnodeX" << std::endl;
        for(int i=0;i<=ppts-1;i++){
            for(int j=0;j<(deg+1);j++){
                std::cout << VnodeX(i,j) << " ";
            }
            std::cout << std::endl << std::flush;
        }
        std::cout << "VnodeY" << std::endl;
        for(int i=0;i<=ppts-1;i++){
            for(int j=0;j<(deg+1);j++){
                std::cout << VnodeY(i,j) << " ";
            }
            std::cout << std::endl << std::flush;
        }*/
        /*std::cout << "V_Node" << std::endl;
        for(int i=0;i<=ppts-1;i++){
            for(int j=0;j<(deg+1)*(deg+1);j++){
                std::cout << V_node_(i,j) << " ";
            }
            std::cout << std::endl << std::flush;
        }*/
    }
}


// compute V*c restricted to node indices
void Node::setVc_node(const Eigen::VectorXd& c)
{
    if(r_ind_ - l_ind_ > 0) {
        Eigen::VectorXd c_seg = c.segment(l_ind_, r_ind_-l_ind_+1);
        Vc_node_ = V_node_.transpose() * c_seg;
    }
    if(PPointsTree_.size() > 1) {
        Eigen::VectorXd c_seg(PPointsTree_.size());
        int k = 0;
        for (std::vector<Point>::iterator it=PPointsTree_.begin(); it!=PPointsTree_.end(); it++, k++) {
            c_seg(k) = c[it->getId()];
        }
        Vc_node_ = V_node_.transpose() * c_seg;
    }
    //std::cout << "Vc_Node" << std::endl;

    //std::cout << Vc_node_ << std::endl << std::flush;


}

void Node::printree(int n)
{
    std::cout << "Node " << n << std::endl;
    for (std::vector<Point>::iterator it=PPointsTree_.begin(); it!=PPointsTree_.end(); it++) {
        std::cout << it->getId() << " ";
    }
    std::cout << std::endl;
    if (tl_child_!=NULL) {
        std::cout << "tl " << n << " ";
        tl_child_->printree(n+1);
    }
    if (tr_child_!=NULL) {
        std::cout << "tr " << n << " ";
        tr_child_->printree(n+1);
    }
    if (bl_child_!=NULL) {
        std::cout << "bl " << n << " ";
        bl_child_->printree(n+1);
    }
    if (br_child_!=NULL) {
        std::cout << "br " << n << " ";
        br_child_->printree(n+1);
    }
}
