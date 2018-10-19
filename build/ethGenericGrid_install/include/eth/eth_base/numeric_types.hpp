/********************************************************************
  created:	2013/07/23
  created:	23:7:2013   15:04
  filename: 	numeric_types.hpp
  author:		Raffael Casagrande
  
  purpose:	Define common numeric types for the whole project.
*********************************************************************/

#ifndef ETH_BASE_NUMERIC_TYPES_HPP
#define ETH_BASE_NUMERIC_TYPES_HPP

// system includes -------------------------------------------------------------
#include <cstddef>
#include <limits>
#include <complex>
#include <boost/limits.hpp>
#include <Eigen/Dense>

namespace eth {
  namespace base {

    //! define an unsigned type
    typedef unsigned unsigned_t;
    //! define an signed type
    typedef int      signed_t;
    //! define a type to use for coordinates:
    typedef double coord_t;

    template< typename T, int ROWS, int COLS >
    using eigen_matrix_t = Eigen::Matrix< T, ROWS, COLS >;

    //! Standard fixed size matrix.
    template<int ROWS, int COLS>
    using fixedMatrix_t = Eigen::Matrix<coord_t,ROWS,COLS>;

    //! Standard fixed size array.
    template<int ROWS, int COLS>
    using fixedArray_t = Eigen::Array<coord_t,ROWS,COLS>;

    //! Standard dynamic size matrix
    template< typename T >
    using dynamicMatrix_t = Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic>;

    //! Standard dynamic size array
    template< typename T >
    using dynamicArray_t = Eigen::Array<T,Eigen::Dynamic,Eigen::Dynamic>;

    //! Standard double dynamic size matrix
    typedef dynamicMatrix_t<double> dynamicMatrixD_t;
    //! Standard complex dynamic size matrix
    typedef dynamicMatrix_t<std::complex<double>> dynamicMatrixC_t;

    /**
     * \brief Query if a coordinate type is equal to zero (absolutely)
     * \param value The value.
     * \return  true if zero, false if not.
     */
    inline bool isZero(coord_t value) {
      return value < 1e-8;
    }

    template<int ROWS, int COLS>
    inline bool isMatrixEqual(const Eigen::Matrix<coord_t,ROWS,COLS>& lhs,
                              const Eigen::Matrix<coord_t,ROWS,COLS>& rhs) 
    {
      if(isZero(lhs.norm()) && isZero(rhs.norm())) return true;
      return isZero((rhs-lhs).norm()/(rhs.norm()+lhs.norm()));
    }

    ////! compile time power
    //template< int N >
    //struct Power
    //{
    //  static const long value = Power< N-1 >::value * 10;
    //};
    //
    //template< >
    //struct Power< 0 >
    //{
    //  static const long value = 1;
    //};

    ////! given an exponent EXP provide an eps = 10^(-|EXP|)
    //template< int EXP >
    //struct Eps {
    //  static const int    EXP_ = ( EXP < 0 ? -EXP : EXP );
    //  static const long   denom;
    //  static const double eps;
    //};

    //template< int EXP >
    //const long Eps< EXP >::denom = Power< EXP >::value;

    //template< int EXP >
    //const double Eps< EXP >::eps = 
    //  1. / static_cast< double >( Eps< EXP >::denom );


    ////! for an eps = 10^-|EXP| check whether a number is "zero"
    //template< int EXP, typename T  >
    //struct IsZero
    //{
    //  bool operator()( const T& number )
    //  {
    //    const double bound = Eps<EXP>::eps;
    //    return ( std::abs( number ) < bound ? true : false );
    //  }
    //};

    ////! convenience function for double numbers
    //template< int EXP >
    //bool is_zero( const double number ) 
    //{ 
    //  return IsZero<EXP,double>()( number ); 
    //}

    ////! convenience function for complex numbers
    //template< int EXP >
    //bool is_zero( const std::complex<double> number ) 
    //{
    //  const double abs_ = std::abs( number );
    //  return IsZero<EXP,double>()( abs_ );
    //}

  } // end namespace base
} // end namespace eth

#endif // ETH_BASE_NUMERIC_TYPES_HPP
