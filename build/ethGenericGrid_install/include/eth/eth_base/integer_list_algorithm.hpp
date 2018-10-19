//------------------------------------------------------------------------------
// <preamble>
// </preamble>
//------------------------------------------------------------------------------

//! @file   integer_list_algorithm.hpp
//! @author Lars Kielhorn
//! @date   2013


#ifndef ETH_GRID_INTEGER_LIST_ALGORITHM_HPP
#define ETH_GRID_INTEGER_LIST_ALGORITHM_HPP

// system includes -------------------------------------------------------------
#include <tuple>

//------------------------------------------------------------------------------
namespace eth {
  namespace base {



    //--------------------------------------------------------------------------
    /** \struct SumIntegerList
     * \brief Compile time summation of an integer list
     */
    struct SumIntegerList
    {
      //! compute \f$ \sum_{i=0}^{N} t_i \,, N = \dim(t)-1 \f$
      template< typename INTEGER_LIST >      
      constexpr static int compute( )
      {
        return Helper_< INTEGER_LIST::size( )-1, INTEGER_LIST >::eval( );
      }
      
    private:
      template< int N, typename IL >
      struct Helper_
      {
        constexpr static int eval( )
        { return IL::template get<N>() + Helper_<N-1,IL>::eval(); }
      };
      
      
      template< typename IL >
      struct Helper_< 0, IL >
      {
        constexpr static int eval( )
        { return IL::template get< 0 >( ); }
      };
      
    };



    //--------------------------------------------------------------------------
    /** \struct DotIntegerList
     * \brief Compile time dot product of two integer lists
     */
    struct DotIntegerList
    {
      //! compute \f$ \sum_{i=0}^{N} t_{1i} t_{2i} \,, N = \dim(t_1)-1 = \dim(t_2)-1 \f$
      template< typename INTEGER_LIST_1,
                typename INTEGER_LIST_2 >
      constexpr static int compute( )
      {
        static_assert( INTEGER_LIST_1::size() == INTEGER_LIST_2::size(),
                       "Integer lists must be of equal length!" );
        return Helper_< INTEGER_LIST_1::size()-1,
                        INTEGER_LIST_1,INTEGER_LIST_2>::eval( );
      }
      
    private:
      template< int N, typename IL_1, typename IL_2 >
      struct Helper_
      {
        constexpr static int eval( )
        {
          return IL_1::template get<N>( ) * IL_2::template get<N>( )
            + Helper_< N-1, IL_1, IL_2 >::eval( );
        }
      };
      
      
      template< typename IL_1, typename IL_2 >
      struct Helper_< 0, IL_1, IL_2 >
      {
        constexpr static int eval( )
        {
          return IL_1::template get<0>() * IL_2:: template get<0>();
        }
      };
      
    };
    

  } // end namespace helper
} // end namespace betl2

#endif // BETL2_UTILS_INTEGER_LIST_ALGORITHM_HPP
