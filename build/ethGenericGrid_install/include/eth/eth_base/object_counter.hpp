//------------------------------------------------------------------------------
// <preamble>
// </preamble>
//------------------------------------------------------------------------------

//! @file   object_counter.hpp
//! @author Lars Kielhorn
//! @date   2013

#ifndef ETH_BASE_ENTITYCOUNTER_HPP
#define ETH_BASE_ENTITYCOUNTER_HPP

// system includes -------------------------------------------------------------
#include <cstddef>
#include <vector>

// own includes ----------------------------------------------------------------
#include "numeric_types.hpp"

//------------------------------------------------------------------------------
namespace eth {
  namespace base {

    template< typename T, typename INDEX_T = eth::base::unsigned_t >
    class EntityCounter;

  } // end namespace base
} // end namespace eth


//------------------------------------------------------------------------------
/** \class EntityCounter
 * \brief The EntityCounter my be used in conjunction with CRTP in order 
 * to automatically keep track of the number of created entitys
 *
 * \tparam T this represents any class which derives from EntityCounter
 * \tparam INDEX_T the index type to be used for counting (default eth::base::unsigned_t)
 */
template< typename T, typename INDEX_T >
class eth::base::EntityCounter {
public:
  typedef std::vector< INDEX_T > vector_t;

private:
  static vector_t       alive_;
  static const unsigned RefElTypesSize_ = 8;

protected:
  EntityCounter( eth::base::RefElType t ) 
  {
    ++( alive_[ static_cast<unsigned>(t) ] );
  }
  
  virtual ~EntityCounter( )
  {
    /* empty */
  }
  
  void destruct_( eth::base::RefElType t )
  {
    --( alive_[ static_cast<unsigned>(t) ] );
  }

public:
  static void resetCounter( )
  {
    for( auto& a : alive_ ) a = 0;
  }

  static void resetCounter( eth::base::RefElType t )
  {
    alive_[ static_cast<unsigned>(t) ] = 0;
  }

  static INDEX_T howMany( eth::base::RefElType t )
  {
    return alive_.at( static_cast<unsigned>(t) );
  }

};

//------------------------------------------------------------------------------
// instantiate static data
template< typename T, typename INDEX_T >
typename eth::base::EntityCounter<T,INDEX_T>::vector_t 
eth::base::EntityCounter<T,INDEX_T>::alive_( eth::base::EntityCounter<T,INDEX_T>::RefElTypesSize_, 0 );

#endif // ETH_BASE_ENTITYCOUNTER_HPP
