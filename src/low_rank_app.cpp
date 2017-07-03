#include "../include/low_rank_app.hpp"
#include "../include/ctree.hpp"
#include "../include/kernel.hpp"

class LowRankApp {

public:

    LowRankApp( Kernel kernel, const Eigen::VectorXd& x, const Eigen::VectorXd& y ):
        kernel_(kernel), Tx_( Tx(x) ), Ty_( Ty(y) )
    { }

    // Approximate matrix-vector multiplication
    Eigen::VectorXd mvProd( const Eigen::VectorXd& c, double eta, unsigned deg )
    {
        // V-matrices
        Tx_.add_V(deg);
        Ty_.add_V(deg);

        // Product V*c restricted to indices of each cluster
        Ty_.vc_mult(c);

        // Admissible partition consisting of near and far field
        Tx_.near_far(eta, Ty_);

        // Compute far-field contribution
        Eigen::VectorXd fapprox = Eigen::VectorXd::Zero(c.size());
        ff_contribution(fapprox, Tx.get_root(), deg);
        // Compute near-field contribution
        nf_contribution(fapprox, Tx.get_root());

        return fapprox;
    }

private:

// far field contribution
void ff_contribution(Eigen::VectorXd& f_, node* x_, int deg,const std::vector<double>& Xvals,const std::vector<double>& Yvals,  double K) {	
  
  if ((*x_).leftchild()!=NULL){ 
    Eigen::VectorXd s=Eigen::VectorXd::Zero(deg+1);	// auxilary variable
    int ixl=(*x_).left_ind(); 	// start index of cluster *x_
    int ixr=(*x_).right_ind(); 	// last index of cluster *x_
    
    // iterate over far field of *x_
    std::vector<node*> ffx=(*x_).get_farf(); // far field of *x_	
    for(std::vector<node*>::iterator iter=ffx.begin(); iter!=ffx.end();++iter){	
      
      int iyl=(**iter).left_ind(); // start index of current cluster in the far field
      int iyr=(**iter).right_ind(); // last index of current cluster in the far field
      BC sigma (Xvals[ixl],Xvals[ixr], Yvals[iyl],Yvals[iyr], deg,K);	 
      const Eigen::MatrixXd& Xx_=sigma.get_matrix(); // the matrix $X_{\sigma,\mu}$
      const Eigen::VectorXd& v_c=(**iter).get_Vc();	
      // V*c restricted to the indices of **iter		
      s=s+Xx_ * v_c; // add the contribution of the block of **iter to s
    };
    
    const Eigen::MatrixXd &Vx_= (*x_).Vv(); 	// $V_{\sigma}$
    f_.segment(ixl,ixr-ixl+1)=f_.segment(ixl,ixr-ixl+1)+Vx_*s; //add the contribution of the far field to f_
    
    node* xl_c=(*x_).leftchild();	// pointer to left child of *x_
    node* xr_c=(*x_).rightchild();	// pointer to right child of *x_
    
    // add contribution of the leaves of *x_
    ff_contribution(f_,xl_c,deg, Xvals,Yvals,K);
    ff_contribution(f_,xr_c,deg, Xvals,Yvals,K);
  };
}

// compute near field contribution
void nf_contribution(Eigen::VectorXd& f_, node* x_,const std::vector<double>& X,const std::vector<double>& Y, const Eigen::VectorXd& c_, double K){
  
  if (x_!=NULL){
    
    int ixl=(*x_).left_ind(); // start index of cluster *x_
    int ixr=(*x_).right_ind(); // last index of cluster *x_
    Kernel G(K);
    
    // iterate over near field of *x_
    std::vector<node*> nfx=(*x_).get_nearf();// near field of *x_
    for(std::vector<node*>::iterator iter=nfx.begin(); iter!=nfx.end();++iter){	
      int iyl=(**iter).left_ind(); // start index of current cluster in the near field
      int iyr=(**iter).right_ind(); // last index of current cluster in the near field
      for(int i=ixl; i<=ixr;++i){	
	for(int j=iyl;j<=iyr;++j){
	  f_(i)=f_(i)+G(X[i],Y[j]) * c_(j); //add near field contribution
	};
      };		
    };
    
    node* xl_c=(*x_).leftchild();// pointer to left child of *x_
    if(xl_c!=NULL){ // if *x_ is not a leaf,...
      node* xr_c=(*x_).rightchild();// pointer to right child of *x_
      // add contribution of the leaves of x_	
      nf_contribution(f_,xl_c,X,Y,c_,K);
      nf_contribution(f_,xr_c,X,Y,c_,K);	
    };
  };
  
}

private:

    Kernel kernel_;
    ctree Tx_;
    ctree Ty_;
};
