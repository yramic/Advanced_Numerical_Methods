//------------------------------------------------------------------------------
// <preamble>
// </preamble>
//------------------------------------------------------------------------------

//! @file   matrix_eth_operations.hpp
//! @author Lars Kielhorn
//! @date   2013

#ifndef ETH_LINALG_MATRIX_ETH_OPERATIONS_HPP
#define ETH_LINALG_MATRIX_ETH_OPERATIONS_HPP

// own includes ----------------------------------------------------------------
#include "is_even.hpp"

//------------------------------------------------------------------------------
namespace eth {
  namespace linalg {
    namespace fixed {

      template< typename SCALAR, 
                std::size_t ROWS, 
                std::size_t DIM, 
                std::size_t COLS,
                enum Storage S >
      Matrix< SCALAR, ROWS, COLS, S > 
      operator*( const Matrix< SCALAR, ROWS, DIM , S >& A,
                 const Matrix< SCALAR, DIM , COLS, S >& B );
    } // end namespace fixed
  } // end namespace linalg
} // end namespace eth


//------------------------------------------------------------------------------
namespace eth {
  namespace linalg {
    namespace fixed {

      namespace /* anonymous */ {

        //----------------------------------------------------------------------
        template< bool COND > struct MatMult { };



        //----------------------------------------------------------------------
        template< > 
        struct MatMult< true > 
        {
          template< typename SCALAR,
                    std::size_t ROWS,
                    std::size_t DIM,
                    std::size_t COLS,
                    enum Storage S >
          Matrix< SCALAR, ROWS, COLS, S >
          operator()( const Matrix< SCALAR, ROWS, DIM , S >& A,
                      const Matrix< SCALAR, DIM , COLS, S >& B ) const
          {
            Matrix< SCALAR, ROWS, COLS, S > C( SCALAR(0) );
            for( std::size_t i = 0; i < ROWS; ++i )
              for( std::size_t j = 0; j < COLS; ++j ) {
                for( std::size_t k = 0; k < DIM ; k = k + 2 ) {
                  C( i, j ) += A( i, k   ) * B( k  , j );
                  C( i, j ) += A( i, k+1 ) * B( k+1, j );
                }
              }
            return C;
          }
        };



        //----------------------------------------------------------------------
        template< >
        struct MatMult< false >
        {
          template< typename SCALAR,
                    std::size_t ROWS,
                    std::size_t DIM,
                    std::size_t COLS,
                    enum Storage S >
          Matrix< SCALAR, ROWS, COLS, S >
          operator()( const Matrix< SCALAR, ROWS, DIM , S >& A,
                      const Matrix< SCALAR, DIM , COLS, S >& B ) const
          {
            Matrix< SCALAR, ROWS, COLS, S > C( SCALAR(0) );
            for( std::size_t i = 0; i < ROWS; ++i )
              for( std::size_t j = 0; j < COLS; ++j ) {
                C( i, j ) += A( i, 0 ) * B( 0, j );
                for( std::size_t k = 1; k < DIM ; k = k + 2 ) {
                  C( i, j ) += A( i, k   ) * B( k  , j );
                  C( i, j ) += A( i, k+1 ) * B( k+1, j );
                }
              }
            return C;
          }
        };

      } // end namespace anonymous
  

      //------------------------------------------------------------------------
      template< typename SCALAR,
                std::size_t ROWS,
                std::size_t DIM,
                std::size_t COLS,
                enum Storage S >
      Matrix< SCALAR, ROWS, COLS, S >
      operator*( const Matrix< SCALAR, ROWS, DIM , S >& A,
                 const Matrix< SCALAR, DIM , COLS, S >& B )
      {
        const bool is_even = isEven< DIM >::value;
        return MatMult< is_even >()( A, B );
      }

    } // end namespace fixed
  } // end namespace linalg
} // end namespace eth


#endif // ETH_LINALG_MATRIX_ETH_OPERATIONS_HPP
