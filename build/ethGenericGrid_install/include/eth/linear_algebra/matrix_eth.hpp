//------------------------------------------------------------------------------
// <preamble>
// </preamble>
//------------------------------------------------------------------------------

//! @file   matrix_eth.hpp
//! @author Lars Kielhorn
//! @date   2013

// system includes -------------------------------------------------------------
#include <array>
#include <cstddef>
#include <cassert>
#include <iostream>

#ifndef ETH_LINALG_ETH_MATRIX_HPP
#define ETH_LINALG_ETH_MATRIX_HPP


//------------------------------------------------------------------------------
namespace eth {
  namespace linalg {
    namespace fixed {
    
      enum class Storage : std::size_t { rowWise, colWise };
      
      template< typename SCALAR, std::size_t ROWS, std::size_t COLS,
                enum Storage S = Storage::colWise >
      class Matrix;
      
    } // end namespace fixed
  } // end namespace linalg
} // end namespace eth



//------------------------------------------------------------------------------
namespace /* anonymous */ {

  //----------------------------------------------------------------------------
  template< enum eth::linalg::fixed::Storage S, 
            std::size_t ROWS, std::size_t COLS > 
  struct StorageTraits { };

    
  //----------------------------------------------------------------------------
  template< std::size_t ROWS, std::size_t COLS >
  struct StorageTraits< eth::linalg::fixed::Storage::colWise, ROWS, COLS > {
    std::size_t operator()( std::size_t i, std::size_t j ) const
    {
      assert( i < ROWS );
      assert( j < COLS );
      return j * ROWS + i;
    }
  };


  //----------------------------------------------------------------------------
  template< std::size_t ROWS, std::size_t COLS >
  struct StorageTraits< eth::linalg::fixed::Storage::rowWise, ROWS, COLS > {
    std::size_t operator()( std::size_t i, std::size_t j ) const
    {
      assert( i < ROWS );
      assert( j < COLS );
      return i * COLS + j;
    }
  };

} // end namespace anonymous

//------------------------------------------------------------------------------
template< typename SCALAR, std::size_t ROWS, std::size_t COLS,
          enum eth::linalg::fixed::Storage S >
class eth::linalg::fixed::Matrix {

  template< typename T, std::size_t R, std::size_t C, enum Storage ST >
  friend
  std::ostream& operator<<( std::ostream& out, 
                            const Matrix< T,R,C,ST >& m );

private:
  typedef Matrix< SCALAR, ROWS, COLS > self_t;
  typedef StorageTraits< S, ROWS, COLS > Traits;

  typedef std::array< SCALAR, ROWS * COLS > array_t;

  array_t data_;

public:
  static const std::size_t NoOfRows = ROWS;
  static const std::size_t NoOfCols = COLS;

public:
  typedef SCALAR        numeric_t;
  typedef SCALAR*       pointer_t;
  typedef const SCALAR* const_pointer_t;
  typedef typename array_t::iterator iterator;
  typedef typename array_t::const_iterator const_iterator;

public:
  Matrix( std::size_t dummy_rows, std::size_t dummy_cols ) 
  { 
    assert( dummy_rows <= ROWS );
    assert( dummy_cols <= COLS );
  }

  Matrix( ) { /* empty */ }

  Matrix( numeric_t val ) { data_.fill( val ); }

  numeric_t& operator()( std::size_t i, std::size_t j ) 
  {
    return data_[ Traits()( i, j ) ];
  }

  numeric_t operator()( std::size_t i, std::size_t j ) const
  {
    return data_[ Traits()( i, j ) ];
  }

  pointer_t       getData( )       { return &data_[ 0 ]; }
  const_pointer_t getData( ) const { return &data_[ 0 ]; }

  iterator       begin( )       { return data_.begin( ); }
  iterator       end  ( )       { return data_.end  ( ); }
  const_iterator begin( ) const { return data_.begin( ); }
  const_iterator end  ( ) const { return data_.end  ( ); }

private:
  std::ostream& write_( std::ostream& out ) const
  {
    out << "[" << ROWS << " x " << COLS << "]:" << std::endl;
    for( std::size_t i = 0; i < ROWS; ++i ) {
      for( std::size_t j = 0; j < COLS; ++j )
        out << data_[ Traits()( i, j ) ] << ( j != COLS-1 ? ", " : " " );
      out << std::endl;
    }
    return out;
  }
};


//------------------------------------------------------------------------------
namespace eth {
  namespace linalg {
    namespace fixed {
      template< typename T, std::size_t R, std::size_t C, enum Storage ST >
      std::ostream& operator<<( std::ostream& out, const Matrix<T,R,C,ST>& m )
      {
        return m.write_( out );
      }
    }
  }
}

//------------------------------------------------------------------------------
// a matrix class without linear algebra operations does not make sense
//  --> include them!!!
#include "matrix_eth_operations.hpp"

#endif // ETH_LINALG_ETH_MATRIX_HPP
