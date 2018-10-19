//------------------------------------------------------------------------------
// <preamble>
// </preamble>
//------------------------------------------------------------------------------

//! @file   integer_list.hpp
//! @author Lars Kielhorn
//! @date   2013

#ifndef ETH_BASE_INTEGER_LIST_HPP
#define ETH_BASE_INTEGER_LIST_HPP


//------------------------------------------------------------------------------
namespace eth {
  namespace base {

    //--------------------------------------------------------------------------
    /** \class integer_list
     * \brief Integer lists are tuple-like structures which are constructed 
     * by an arbitrary number of integers as non-type template arguments. 
     * Simple compile-time algorithms on basis of integer lists can be found in 
     * integer_list_algorithm.hpp.
     * 
     */
    template< int HEAD, int... TAIL >
    class integer_list {
    public:
      /// get N-th element
      template< int N >
      constexpr static int get( )
      {
        static_assert( N < size( ), 
                       "Given value exceeds the tuple's size!" );
        return get_<N,HEAD,TAIL...>::get( );    
      }
      /// get the integer list's size
      constexpr static int size( )
      {
        return sizeof...( TAIL ) + 1;
      }

      template< int M >
      constexpr static bool isMember( )
      {
        return isMember_< M, size()-1 >::eval( );
      }


    private:
      template< int N, int H, int... T >
      struct get_
      {
        constexpr static int get( )
        {
          return get_<N-1,T...>::get( );
        }
      };
  
      template< int H, int... T >
      struct get_<0,H,T...>
      {
        constexpr static int get( )
        {
          return H;
        }
      };

      //------------------------------------------------------------------------
      template< int M, int N >
      struct isMember_
      { 
        constexpr static bool eval( )
        {
          return ( M == get<N>() ? true : isMember_<M,N-1>::eval() );
        }
      };

      template< int M >
      struct isMember_< M, 0 >
      {
        constexpr static bool eval( )
        {
          return ( M == get<0>() ? true : false );
        }
      };

    }; 


    //--------------------------------------------------------------------------
    /** \class integer_list
     * \brief Terminating implementation of the integer_list class
     */
    template< int HEAD >
    class integer_list< HEAD >
    {
    public:
      template< int N >
      constexpr static int get( )
      {
        static_assert( N < 1, 
                       "Given value exceeds the lists's size!" );
        return HEAD;
      }

      constexpr static int size( )
      {
        return 1;
      }

      template< int M >
      constexpr static bool isMember( )
      {
        return ( M == HEAD ? true : false );
      }

    };

  } // end namespace base
} // end namespace eth



#endif // ETH_BASE_INTEGER_LIST_HPP
