/********************************************************************
  created:	2013/06/17
  created:	17:6:2013   17:37
  filename: 	entityIterator.hpp
  author:		Raffael Casagrande
  
  purpose:	Provide an abstract interface for an entity iterator
*********************************************************************/

#ifndef entityIterator_HPP__
#define entityIterator_HPP__

namespace eth{
  namespace grid{
    template<class GRID_TRAITS, class EI_IMPL>
    class EntityIterator;
  }
}

// Implementation
//////////////////////////////////////////////////////////////////////////
#include "entity_pointer.hpp"

/**
 * \brief	Provides an abstract interface to iterate through a collection of 
 * 			Entities (e.g. LeafView, or Hierarchic Iterator etc...)
 * \tparam GRID_TRAITS The GridTraits class (see \ref GridTraitsDoc)
 * \tparam EI_IMPL Underlying EntityIterator implementation class. (Can be 
 * 				   different for codimension and leaf/level/hierarchic iterators)
 * 
 * An EntityIterator derives from a EntityPointer and can always be
 * assigned/compared to an EntityPointer (see below for more details). In 
 * addition to the EntityPointer the EntityIterator provides the operator++
 * (prefix) to move it to the next entity in the set.
 * 
 * As an EntityPointer can only point to entities of a fixed codimension, so does
 * the EntityIterator. The codimension is given be EI_IMPL::codim.
 * 
 * <H3> Relation to EntityPointers </H3>
 * In order for EntityIterators to be comparable to other types of EntityIterators
 * a rather sophisticated inheritance hierarchy must be implemented on
 * both the interface and implementation side. See 
 * eth::grid::EntityPointer  for more details.
 * 
 * <H3>Special Implementation specific functions </H3>
 * The following functions must be implemented solely by the engine class:
 * 
 * Implementation function                                   | Corresponding Interface function(s)
 * ----------------------------------------------------------|-------------------------------------------------
 * void increment()				                     							 | EntityIterator operator++()
 * 
 * The functions which are listed on the rhs of the above table must not be implemented in the engine class.
 */
template<class GRID_TRAITS, class EI_IMPL>
class eth::grid::EntityIterator : public EntityPointer<GRID_TRAITS,EI_IMPL> {
public:
  typedef EntityPointer<GRID_TRAITS,EI_IMPL> base_t;

  /// Construct EntityIterator from underlying implementation (engine pattern)
  explicit EntityIterator(const EI_IMPL& impl) : base_t(impl) {}

  /**
   * \brief	Move the iterator to the next position (Forward iterator, prefix increment)
   * \return	Iterator at new position (++i, prefix increment)
   * \remarks Do not implement this method in engine class, instead provide
   * 			the method void increment();
   */
  EntityIterator& operator++() 
  {
    base_t::impl_.increment();
    return *this;
  }
};

#endif // entityIterator_HPP__
