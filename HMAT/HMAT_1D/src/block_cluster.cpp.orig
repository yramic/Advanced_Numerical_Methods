#include "../include/block_cluster.hpp"
#include "../include/cheby.hpp"


// constructor
BlockCluster::BlockCluster(double xl, double xr, double yl, double yr, unsigned deg, Kernel G):
    xl_(xl), xr_(xr), yl_(yl), yr_(yr), deg_(deg), G_(G), X_(Eigen::MatrixXd::Zero(deg+1,deg+1))
{
    setMatrix();
}


// compute matrix $X_{\sigma,\mu}$
void BlockCluster::setMatrix()
{
    Cheby C_ptx(xl_, xr_, deg_);
    Cheby C_pty(yl_, yr_, deg_);

    Eigen::VectorXd cheb_ptx = C_ptx.getNodes(); // Chebyshew nodes in interval [xl,xr]
    Eigen::VectorXd cheb_pty = C_pty.getNodes(); // Chebyshew nodes in interval [yl,yr]

    for(unsigned i=0; i<=deg_; ++i)
        for(unsigned j=0; j<=deg_; ++j)
            X_(i,j) = G_(cheb_ptx[i],cheb_pty[j]);
}
