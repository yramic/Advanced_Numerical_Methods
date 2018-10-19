//------------------------------------------------------------------------------
// <preamble>
// </preamble>
//------------------------------------------------------------------------------

//! @file   matrix_eigen3.hpp
//! @author Lars Kielhorn
//! @date   2013

#ifndef ETH_LINALG_MATRIX_EIGEN3_HPP
#define ETH_LINALG_MATRIX_EIGEN3_HPP

// eigen3 includes -------------------------------------------------------------
#include <Eigen/Dense>

//------------------------------------------------------------------------------
namespace eth {
  namespace linalg {
    namespace fixed {

      template< typename SCALAR, std::size_t ROWS, std::size_t COLS >
      class Matrix;

      template< typename SCALAR, std::size_t ROWS >
      class Vector;

    }
  }
}

//------------------------------------------------------------------
/** \class Matrix
 * \brief Matrix is as wrapper for the Eigen3-fixed size matrix format
 * with compile time bounds and a column-wise storage scheme
 * 
 * \tparam SCALAR the basic numeric type
 * \tparam ROWS the number of rows
 * \tparam COLS the number of columns
 */
template< typename SCALAR, std::size_t ROWS, std::size_t COLS >
class eth::linalg::fixed::Matrix : public ::Eigen::Matrix< SCALAR, ROWS, COLS >
{
private:
  /** @name derived and base class types */
  //@{
  typedef Matrix< SCALAR, ROWS, COLS >          self_t;
  typedef ::Eigen::Matrix< SCALAR, ROWS, COLS > base_t;
  //@}

public:
  typedef SCALAR numeric_t;
  /** @name static data */
  //@{
  static const std::size_t rows = ROWS; //!< number of rows
  static const std::size_t cols = COLS; //!< number of columns
  //@}

public:
  /** @name constructors, assignements, and copies */
  //@{
  //! delegate construction to base class
  Matrix( std::size_t r, std::size_t c ) : base_t( r, c )
  { /* empty */ }
  //! default cnstructor
  Matrix( ) : base_t( ) { /* empty */ }
  //! construct from any matrix expression
  template< typename OTHER_EXPRESSION >
  Matrix( const OTHER_EXPRESSION& cp ) : base_t( cp ) { /* empty */ }
  //! assign with any matrix expression
  template< typename OTHER_EXPRESSION >
  self_t& operator=( const OTHER_EXPRESSION& expr )
  {
    base_t::operator=( expr );
    return *this;
  }
  //@}
};



//------------------------------------------------------------------
/** \class Vector
 * \brief Vector is as wrapper for the Eigen3-fixed size matrix format
 * with compile time bounds and a column-wise storage scheme
 * 
 * \tparam SCALAR the basic numeric type
 * \tparam ROWS the number of rows
 * \tparam COLS the number of columns
 */
template< typename SCALAR, std::size_t ROWS >
class eth::linalg::fixed::Vector : public ::Eigen::Matrix< SCALAR, ROWS, 1 >
{
private:
  /** @name derived and base class types */
  //@{
  typedef Vector< SCALAR, ROWS >                self_t;
  typedef ::Eigen::Matrix< SCALAR, ROWS, 1 > base_t;
  //@}

public:
  typedef SCALAR numeric_t;
  /** @name static data */
  //@{
  static const std::size_t rows = ROWS; //!< number of rows
  static const std::size_t cols = 1;    //!< number of columns (=1)
  //@}

public:
  /** @name constructors, assignements, and copies */
  //@{
  //! delegate construction to base class
  Vector( std::size_t r ) : base_t( r ) { /* empty */ }
  //! default cnstructor
  Vector( ) : base_t( ) { /* empty */ }
  //! construct from any matrix expression
  template< typename OTHER_EXPRESSION >
  Vector( const OTHER_EXPRESSION& cp ) : base_t( cp ) { /* empty */ }
  //! assign with any matrix expression
  template< typename OTHER_EXPRESSION >
  self_t& operator=( const OTHER_EXPRESSION& expr )
  {
    base_t::operator=( expr );
    return *this;
  }
  //@}
};


#endif // ETH_LINALG_MATRIX_EIGEN3_HPP
