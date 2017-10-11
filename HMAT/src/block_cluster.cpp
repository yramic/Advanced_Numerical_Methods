#include "../include/block_cluster.hpp"
#include "../include/cheby.hpp"


// constructor for BEM(4d)
BlockCluster::BlockCluster(double x1l, double x1r, double y1l, double y1r, double x2l, double x2r, double y2l, double y2r, unsigned deg, Kernel G):
    x1l_(x1l), x1r_(x1r), y1l_(y1l), y1r_(y1r), x2l_(x2l), x2r_(x2r), y2l_(y2l), y2r_(y2r), deg_(deg), G_(G), X_(Eigen::MatrixXd::Zero((deg+1)*(deg+1),(deg+1)*(deg+1)))
{
    setMatrix();
}

// constructor
BlockCluster::BlockCluster(double xl, double xr, double yl, double yr, unsigned deg, Kernel G):
    xl_(xl), xr_(xr), yl_(yl), yr_(yr), deg_(deg), G_(G), X_(Eigen::MatrixXd::Zero(deg+1,deg+1))
{
    setMatrix();
}

// compute matrix $X_{\sigma,\mu}$
void BlockCluster::setMatrix()
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


    for(unsigned i=0; i<=deg_; ++i)
        for(unsigned j=0; j<=deg_; ++j)
            for(unsigned k=0; k<=deg_; ++k)
                for(unsigned l=0; l<=deg_; ++l)
                    X_(j*(deg_+1)+i,l*(deg_+1)+k) = G_(cheb_pt1x[i],cheb_pt1y[j],cheb_pt2x[k],cheb_pt2y[l]);

}


// 1 time for 2 boxes m(i,j)=
