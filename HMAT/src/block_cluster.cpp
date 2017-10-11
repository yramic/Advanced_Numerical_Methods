#include "../include/block_cluster.hpp"
#include "../include/cheby.hpp"


// constructor
BlockCluster::BlockCluster(double xl, double xr, double yl, double yr, unsigned deg, Kernel G):
    xl_(xl), xr_(xr), yl_(yl), yr_(yr), deg_(deg), G_(G), X_(Eigen::MatrixXd::Zero((deg+1)*(deg+1),(deg+1)*(deg+1)))
{
    setMatrix();
}


// compute matrix $X_{\sigma,\mu}$
void BlockCluster::setMatrix()
{
    // first box
    //Cheby C_pt1x(x1l_, x1r_, deg_);
    //Cheby C_pt1y(y1l_, y1r_, deg_);
    // second box
    //Cheby C_pt2x(x2l_, x2r_, deg_);
    //Cheby C_pt2y(y2l_, y2r_, deg_);

    //Eigen::VectorXd cheb_pt1x = C_pt1x.getNodes(); // Chebyshew nodes in interval [xl,xr]
    //Eigen::VectorXd cheb_pt1y = C_pt1y.getNodes(); // Chebyshew nodes in interval [yl,yr]

    //Eigen::VectorXd cheb_pt2x = C_pt2x.getNodes(); // Chebyshew nodes in interval [xl,xr]
    //Eigen::VectorXd cheb_pt2y = C_pt2y.getNodes(); // Chebyshew nodes in interval [yl,yr]


    /*for(unsigned i=0; i<=deg_; ++i)
        for(unsigned j=0; j<=deg_; ++j)
            for(unsigned k=0; k<=deg_; ++k)
                for(unsigned l=0; l<=deg_; ++l)
                    X_(j*(deg+1)+i,l*(deg+1)+k) = G_(cheb_pt1x[i],cheb_pt1y[j],cheb_pt2x[k],cheb_pt2y[l]);
    */
}


// 1 time for 2 boxes m(i,j)=
