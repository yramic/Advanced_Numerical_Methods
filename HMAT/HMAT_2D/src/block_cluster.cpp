#include "../include/block_cluster.hpp"
#include "../include/cheby.hpp"
#include <iostream>

// constructor for BEM(4d)
/*BlockCluster::BlockCluster(double x1l, double x1r, double y1l, double y1r, double x2l, double x2r, double y2l, double y2r, unsigned deg, Kernel4D* G):
    x1l_(x1l), x1r_(x1r), y1l_(y1l), y1r_(y1r), x2l_(x2l), x2r_(x2r), y2l_(y2l), y2r_(y2r), deg_(deg), G4_(G), X_(Eigen::MatrixXd::Zero((deg+1)*(deg+1),(deg+1)*(deg+1)))
{
    setMatrix4D(G);
}

BlockCluster::BlockCluster(double x1l, double x1r, double y1l, double y1r, double x2l, double x2r, double y2l, double y2r, unsigned deg, PolynomialKernel* G):
    x1l_(x1l), x1r_(x1r), y1l_(y1l), y1r_(y1r), x2l_(x2l), x2r_(x2r), y2l_(y2l), y2r_(y2r), deg_(deg), GP_(G), X_(Eigen::MatrixXd::Zero((deg+1)*(deg+1),(deg+1)*(deg+1)))
{
    setMatrix4D(G);
}
BlockCluster::BlockCluster(double x1l, double x1r, double y1l, double y1r, double x2l, double x2r, double y2l, double y2r, unsigned deg, ConstantKernel* G):
    x1l_(x1l), x1r_(x1r), y1l_(y1l), y1r_(y1r), x2l_(x2l), x2r_(x2r), y2l_(y2l), y2r_(y2r), deg_(deg), GC_(G), X_(Eigen::MatrixXd::Zero((deg+1)*(deg+1),(deg+1)*(deg+1)))
{
    setMatrix4D(G);
}*/
BlockCluster::BlockCluster(double x1l, double x1r, double y1l, double y1r, double x2l, double x2r, double y2l, double y2r, unsigned deg, Kernel* G):
    x1l_(x1l), x1r_(x1r), y1l_(y1l), y1r_(y1r), x2l_(x2l), x2r_(x2r), y2l_(y2l), y2r_(y2r), deg_(deg), X_(Eigen::MatrixXd::Zero((deg+1)*(deg+1),(deg+1)*(deg+1)))
{
    setMatrix4D(G);
}
// constructor
BlockCluster::BlockCluster(double xl, double xr, double yl, double yr, unsigned deg, Kernel2D G):
    xl_(xl), xr_(xr), yl_(yl), yr_(yr), deg_(deg), G2_(G), X_(Eigen::MatrixXd::Zero(deg+1,deg+1))
{
    setMatrix2D();
}

// compute matrix $X_{\sigma,\mu}$
void BlockCluster::setMatrix4D(Kernel* G)
{
    // first box
    Cheby C_pt1x(x1l_, x1r_, deg_);
    Cheby C_pt1y(y1l_, y1r_, deg_);
    // second box
    Cheby C_pt2x(x2l_, x2r_, deg_);
    Cheby C_pt2y(y2l_, y2r_, deg_);

    Eigen::VectorXd cheb_pt1x = C_pt1x.getNodes(); // Chebyshew nodes in interval [x1l,x1r]
    Eigen::VectorXd cheb_pt1y = C_pt1y.getNodes(); // Chebyshew nodes in interval [y1l,y1r]

    Eigen::VectorXd cheb_pt2x = C_pt2x.getNodes(); // Chebyshew nodes in interval [x2l,x2r]
    Eigen::VectorXd cheb_pt2y = C_pt2y.getNodes(); // Chebyshew nodes in interval [y2l,y2r]

    for(int i=0; i<=deg_; i++){
        for(int j=0; j<=deg_; j++){
            for(int k=0; k<=deg_; k++){
                for(int l=0; l<=deg_; l++){
                    X_(i*(deg_+1)+j,k*(deg_+1)+l) = (*G)(cheb_pt1x[i],cheb_pt1y[j],cheb_pt2x[k],cheb_pt2y[l]);
                }
            }
        }
    }
}

void BlockCluster::setMatrix2D()
{
    Cheby C_ptx(xl_, xr_, deg_);
    Cheby C_pty(yl_, yr_, deg_);

    Eigen::VectorXd cheb_ptx = C_ptx.getNodes(); // Chebyshew nodes in interval [xl,xr]
    Eigen::VectorXd cheb_pty = C_pty.getNodes(); // Chebyshew nodes in interval [yl,yr]

    for(unsigned i=0; i<=deg_; ++i)
        for(unsigned j=0; j<=deg_; ++j)
            X_(i,j) = G2_(cheb_ptx[i],cheb_pty[j]);
}


// 1 time for 2 boxes m(i,j)=
