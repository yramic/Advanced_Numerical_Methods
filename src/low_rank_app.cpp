#include "../include/low_rank_app.hpp"
#include "../include/BC.hpp"
#include "../include/ctree.hpp"
#include "../include/kernel.hpp"
#include "../include/node.hpp"


LowRankApp::LowRankApp( Kernel kernel, const Eigen::VectorXd& x, const Eigen::VectorXd& y ):
    kernel_(kernel), Tx_( Tx(x) ), Ty_( Ty(y) )
{ }


// Approximate matrix-vector multiplication
Eigen::VectorXd LowRankApp::mvProd( const Eigen::VectorXd& c, double eta, unsigned deg )
{
    // V-matrices
    Tx_.add_V(deg);
    Ty_.add_V(deg);

    // Product V*c restricted to indices of each cluster
    Ty_.vc_mult(c);

    // Admissible partition consisting of near and far field
    Tx_.near_far(eta, Ty_);

    // Compute far-field contribution
    Eigen::VectorXd f_approx = Eigen::VectorXd::Zero(c.size());
    ff_contribution(f_approx, Tx.get_root(), deg);
    // Compute near-field contribution
    nf_contribution(f_approx, Tx.get_root(), c);

    return f_approx;
}


// Compute far-field contribution
void LowRankApp::ff_contribution( Eigen::VectorXd& f, node* tx, int deg )
{
    if((*tx).leftchild() != NULL) {

        Eigen::VectorXd s = Eigen::VectorXd::Zero(deg+1); // auxilary variable
        int ixl = (*tx).left_ind(); // start index of cluster *tx
        int ixr = (*tx).right_ind(); // last index of cluster *tx

        // Iterate over far field of *tx_root
        std::vector<node*> ffx = (*tx).get_farf(); // far field of *tx_root
        for(std::vector<node*>::iterator iter=ffx.begin(); iter!=ffx.end(); ++iter) {

            int iyl = (**iter).left_ind(); // start index of current cluster in the far field
            int iyr = (**iter).right_ind(); // last index of current cluster in the far field
            BC sigma(Tx_.getVals()[ixl], Tx_.getVals()[ixr], Ty_.getVals()[iyl], Ty_.getVals()[iyr], deg, kernel_);
            Eigen::MatrixXd Xx_ = sigma.get_matrix(); // matrix $X_{\sigma,\mu}$
            Eigen::VectorXd Vc = (**iter).get_Vc();
            // V*c restricted to the indices of **iter
            s += Xx * Vc; // add contribution of block **iter to s
        };

        Eigen::MatrixXd Vx = (*tx).Vv(); // $V_{\sigma}$
        f.segment(ixl, ixr-ixl+1) = f.segment(ixl, ixr-ixl+1) + Vx*s; // add contribution of far field to f_

        node* xl_c = (*tx).leftchild();  // pointer to left  child of *tx
        node* xr_c = (*tx).rightchild(); // pointer to right child of *tx

        // add contribution of leaves of *tx
        ff_contribution(f, xl_c, deg);
        ff_contribution(f, xr_c, deg);
    }
}


// Compute near-field contribution
void LowRankApp::nf_contribution( Eigen::VectorXd& f, node* tx, const Eigen::VectorXd& c_ )
{
    if(Tx != NULL) {

        unsigned ixl=(*tx).left_ind(); // start index of cluster *tx
        unsigned ixr=(*tx).right_ind(); // last index of cluster *tx

        // Iterate over near field of *tx
        std::vector<node*> nfx = (*tx).get_nearf(); // near field of *tx
        for(std::vector<node*>::iterator iter=nfx.begin(); iter!=nfx.end(); ++iter) {
            unsigned iyl = (**iter).left_ind(); // start index of current cluster in the near field
            unsigned iyr = (**iter).right_ind(); // last index of current cluster in the near field
            for(unsigned i=ixl; i<=ixr; ++i) {
                for(unsigned j=iyl; j<=iyr; ++j) {
                    f(i) += kernel_(Tx_.getVals()[i], Ty_.getVals()[j]) * c_(j); // add near field contribution
                }
            }
        }

        node* xl_c = (*tx).leftchild(); // pointer to left child of *tx
        if(xl_c != NULL) { // if *tx is not a leaf
            node* xr_c = (*tx).rightchild(); // pointer to right child of *tx

            // add contribution of leaves of *tx
            nf_contribution(f, xl_c, c_);
            nf_contribution(f, xr_c, c_);
        }
    }
}
