#ifndef BEM_TRANSMISSION_SYSTEM_MATRIX_HPP
#define BEM_TRANSMISSION_SYSTEM_MATRIX_HPP 

// Eigen includs ---------------------------------------------------------------
#include <Eigen/Dense>
#include <Eigen/Sparse>

//----------------------------------------------------------------------------
class TransmissionSystemMatrix 
{
public:
  typedef Eigen::MatrixXd  matrix_t;
  typedef Eigen::VectorXd  vector_t;
    
private:
  const matrix_t& W_;
  const matrix_t& K_;
  const matrix_t& V_;
  const double & alpha_;
  matrix_t A_;

public:
  explicit TransmissionSystemMatrix( const matrix_t&  V, const matrix_t&  K,
			 const matrix_t&  W, const double & alpha)
    : V_( V )
    , K_( K )
    , W_( W ),
      alpha_(alpha)
  {
    // Check non-zero alpha
    if( fabs(alpha) < 1e-12 ) {
      std::cerr <<  "alpha cannot be zero!" << std::endl;
      exit( -1 );
    }
    
    this -> check_matrix_consistencies_();
		
   /// set system matrix to 0. Optional initialization by executing compute()
    A_.resize(0,0);
    A_.setZero();
  }

  /// forbid copies
  TransmissionSystemMatrix( const TransmissionSystemMatrix& ) = delete;

  /// forbid assignments 
  TransmissionSystemMatrix& operator=( const TransmissionSystemMatrix& ) = delete;

  /// destructor
  ~TransmissionSystemMatrix() { /*empty*/ }

    
  void compute( ){
    //construct system matrix using its components (according to 1.13.3)
    A_.conservativeResize(K_.rows() + W_.rows(), K_.cols() + V_.cols());
    A_.block(0, 0, K_.rows(), K_.cols()) = 2.* K_;
    A_.block(0, K_.cols(), V_.rows(), V_.cols()) = -(1./alpha_ + 1)*V_;
    A_.block(K_.rows(), 0, W_.rows(), W_.cols()) = -(alpha_ + 1.)*W_;
    A_.block(K_.rows(), K_.cols(), K_.cols(), K_.rows()) = - 2.*K_.transpose();
  }
	
  const matrix_t& matrix( ) const {
    if ( A_.cols()==A_.rows() && A_.rows()==1 && A_(0,0) == 0. ){
      std::cerr << "Oops, TransmissionSystemMatrix has not yet been computed! "
		<< "Run compute() first!" << std::endl;
      exit(-1);
    }
    return A_;
  }
private:
  void check_matrix_consistencies_( ) const
  {
    if( K_.rows() != V_.rows() ) { std::cerr << "rows(K) != rows(V). Abort!"
					     << std::endl; exit(-1); }
    if( K_.transpose().rows() != W_.rows() ){ std::cerr << "cols(K) != rows(W). Abort!"
							 << std::endl; exit(-1); }
    if( K_.cols() != W_.cols() )        { std::cerr << "cols(K) != cols(W). Abort!"
						    << std::endl; exit(-1); }
    if( K_.transpose().cols() != V_.cols() ){ std::cerr << "rows(K) != cols(V). Abort!"
							<< std::endl; exit(-1); }

  }


}; // end class TransmissionTransmissionSystemMatrix


#endif // BEM_TRANSMISSION_SYSTEM_MATRIX_HPP
