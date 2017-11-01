#include "../include/block_cluster.hpp"
#include "../include/cheby.hpp"
#include <iostream>

// constructor for BEM(2D)
BlockCluster::BlockCluster(Eigen::VectorXd tk1x, Eigen::VectorXd tk1y, Eigen::VectorXd tk2x, Eigen::VectorXd tk2y, unsigned deg, Kernel* G):
    tk1x_(tk1x), tk1y_(tk1y), tk2x_(tk2x), tk2y_(tk2y), deg_(deg), X_(Eigen::MatrixXd::Zero((deg+1)*(deg+1),(deg+1)*(deg+1)))
{
    setMatrix2D(G);
}

// compute matrix $X_{\sigma,\mu}$
void BlockCluster::setMatrix2D(Kernel* G)
{
    for(int i=0; i<=deg_; i++){
        for(int j=0; j<=deg_; j++){
            for(int k=0; k<=deg_; k++){
                for(int l=0; l<=deg_; l++){
                    X_(i*(deg_+1)+j,k*(deg_+1)+l) = (*G)(tk1x_[i],tk1y_[j],tk2x_[k],tk2y_[l]);
                }
            }
        }
    }
}
