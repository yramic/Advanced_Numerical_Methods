//------------------------------------------------------------------------------
// <preamble>
// </preamble>
//------------------------------------------------------------------------------

//! @file   enum_class_hash.hpp
//! @author Lars Kielhorn
//! @date   2013

#ifndef ETH_BASE_HASH_FUNCTIONS_HPP
#define ETH_BASE_HASH_FUNCTIONS_HPP

// system includes -------------------------------------------------------------
#include <functional>

namespace eth {
  namespace base {

    //------------------------------------------------------------------------
    /** \struct EnumClassHash 
     * \brief This struct may be used in conjunction with unordered 
     * associative containers where the Key is a scoped enumerator.
     *
     * \tparam E the scoped enumerator
     * \tparam T the enumerator's base type (default is int)
     */
    //------------------------------------------------------------------------
    template< typename E, typename T = int >
    struct EnumClassHash {
      
      typename std::hash<T>::result_type operator()( E e ) const {
        return std::hash<T>()( static_cast<T>( e ) );
      }
      
    };

    //--------------------------------------------------------------------------
    /** \struct PairHash
     * \brief If pairs are used as keys for unordered associative containers
     * PairHash may be used as a valid hash function. Implementation details
     * are taken from boost::hash_combine
     */
    //--------------------------------------------------------------------------
    struct PairHash {

      /** the hash-function for pairs
       * \tparam T1 the pair's first argument type, std::hash<T1> must be valid
       * \tparam T2 the pair's second argument type, std::hash<T1> must be valid
       */
      template< typename T1, typename T2 >
      typename std::hash<int>::result_type 
      operator()( const std::pair< T1, T2 >& p ) const 
      {
        typename std::hash<int>::result_type seed = 0;
        const int p1 = static_cast< int >( p.first  );
        const int p2 = static_cast< int >( p.second );
        seed ^= std::hash<int>()( p1 ) + 0x9e3779b9 + (seed<<6) + (seed>>2);
        seed ^= std::hash<int>()( p2 ) + 0x9e3779b9 + (seed<<6) + (seed>>2);
        return seed;
      }

    };


    //--------------------------------------------------------------------------
    /** \struct ContainerPointerHash \brief In rare cases, one might want to
     * use small containers made of pointers as keys for unordered associative
     * containers. ContainerPointerHash accomplishes this task. For this, the
     * container MUST be passed in sorted form.
     * Implementation details are taken from boost::hash_combine
     */
    //--------------------------------------------------------------------------
    struct ContainerPointerHash {
      
      template< typename CONTAINER_T >
      typename std::hash<typename CONTAINER_T::value_type>::result_type
      operator()( const CONTAINER_T& c ) const 
      {
        typedef typename CONTAINER_T::value_type _T;
        // initialize seed
        typename std::hash<_T>::result_type seed = 0;
        // go through container
        for( auto elem : c ) 
          seed ^= std::hash<_T>()(elem) + 0x9e3779b9 + (seed<<6) + (seed>>2);
        return seed;
      }
    };

  } // end namespace base
} // end namespace eth

#endif // ETH_BASE_HASH_FUNCTIONS_HPP
