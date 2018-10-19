/********************************************************************
  created:	2013/06/18
  created:	18:6:2013   15:13
  filename: 	entityPointer.hpp
  author:		Raffael Casagrande
  
  purpose:	Abstract interface for an EntityPointer
*********************************************************************/


#ifndef entityPointer_HPP__
#define entityPointer_HPP__

namespace eth{
  namespace grid{
    template<class GRID_TRAITS, class EP_IMPL>
    class EntityPointer;
  }
}


// Implementation
//////////////////////////////////////////////////////////////////////////

#include "entity.hpp"

/**
 * \brief	Interface for an Entity Pointer.
 * \tparam EP_IMPL The underlying implementation class. EP_IMPL::codimen
 * 				   specifies the codimension to which this EntityPointer is
 * 				   pointing.
 * \tparam GRID_TRAITS The GridTraits class (see \ref GridTraitsDoc )
 * 
 * An EntityPointer points typically to an Entity object. Because Entities are
 * templatized on their Codimension, one Entity Pointer can only point to
 * entities of a given codim (determined through the EP_IMPL::codimen variable).
 * \remarks EntityPointers are always immutable, that is the referenced Entity
 * 			can never be changed.
 * \remarks This interface uses the Engine pattern.
 * \remarks In contrast to an Entity object, a EntityPointer supports copying and
 * 			assignment operations.
 * 
 * <H3>Relation to EntityIterators </H3>
 * EntityIterators derive from EntityPointers and add the increment operator++.
 * This inheritance hierarchy should also be reflected in the actual implementation
 * classes which are wrapped by EntityIterator and EntityPointer
 * classes. The following image shows how EntityPointners and EntityIterators
 * are connected on the implementation and interface side:
 * \image html EntityPointerIterator.png "Relation of EntityPointers and EntityIterators"
 * 
 * Note that,
 * - An EntityPointer can be initialized with a EntityIterator, but an EntityIterator cannot be initialized by a EntityPointer.   
 * - An EntityIterator can be compared with any other EntityIterator/EntityPointer. 
 * - All methods in the EntityPointerImplementation class are NON-VIRTUAL!
 * - All EntityIterator implementations must inherit from the same EntityPointerImplementation (per codim).  --> There is only one EntityPointer implementation per codimension but multiple EntityIteratorImplementations per codimension.
 * - The type EP_IMPL2 which appears in method of the EntityPointer interface must inherit from EntityPointerImplementation.
 
 *
 * The advantage of this approach is, that it is possible to assign and compare
 * e.g. a  * LevelIterator to an EntityPointer. On the other hand the implementation  
 * must use the same EntityPointer classes for all EntityIterator implementations,
 * which may not seem to make sense at first sight. However it is much more
 * user-friendly in the sense that any EntityPointer can be compared to any other
 * EntityPointer.
 * 
 * <H3> Implementation specific functions </H3>
 * For clarity we list here all the functions which must be implemented by the
 * engine class. Due to the above inheritance hierarchy, the signatures are quite
 * a bit different from the signatures of the interface.
 * 
 * Implementation function                                      | Corresponding Interface function(s)
 * -------------------------------------------------------------|-------------------------------------------------
 * EP_IMPL(const EP_IMPL& other)  (constructor)                 | template<class EP_IMPL2> EntityPointer(const EntityPointer<gridTraits_t,EP_IMPL2>& other)
 * void assign(const EP_IMPL& rhs)							              	| template<class EP_IMPL2> EntityPointer operator=(const EntityPointer<gridTraits_t,EP_IMPL2>& rhs)
 * const gridTraits_t::Entity<codim>& dereference() const	     	| operator* and operator-> (dereferencing)
 * bool operator==(const EP_IMPL& rhs) const				          	| operator== and operator!= (comparison)
 * int level() const											                      | int level() const
 * 
 * where EP_IMPL means the concrete engine class. Note that the copy constructor,
 * assign and the operator== method in the engine class are NOT templated.
 */
template<typename GRID_TRAITS, typename EP_IMPL>
class eth::grid::EntityPointer {
  // Declare instances of Entity Pointer with arbitrary template arguments as friend.
  template<class,class> friend class EntityPointer;
public:
  /// the traits class which describes the generic Grid Implementation
  typedef GRID_TRAITS gridTraits_t;
  // Underlying implementation class of pointer.
  typedef EP_IMPL impl_t;
  // The codimension of the entities to which the pointer is pointing.
  static const int codim = impl_t::codim;

  typedef typename gridTraits_t::size_type size_type;

  /**
   * \brief	Generic copy constructor: An entity pointer can be initialized
   * 			from a another EntityPointer but also from a LevelIterator or
   * 			a LeafIterator.
   * \tparam the underlying implementation of the other pointer.
   * \param	other	The other EntityPointer.
   * \remarks The underlying implementation must have a constructor which
   * 			takes the other implementation...
   */
  template<class EP_IMPL2>
  EntityPointer(const EntityPointer<gridTraits_t,EP_IMPL2>& other) 
    : impl_(other.impl_) 
  {
    static_assert(std::is_base_of<typename gridTraits_t::template entityPointer_t<codim>,EP_IMPL2>::value,
      "template type EP_IMPL2 must inherit from gridTraits_t::entityPointer_t<codimension>");
  }


  /**
   * \brief	Create pointer from entity.
   * \param	e	The entity to which the pointer shall point.
   * \remarks The EntityPointer implementation must have a constructor taking
   * 			a GRID_TRAITS::entity_t<codim> object.
   */
  // Disable this method, is probably not being used at all.
  //EntityPointer(const Entity<gridTraits_t,codim>& e) : impl_(e.impl_) {}


  /** \brief Wrap an entity pointer around the given implementation
   *  \param impl the implementation to be wrapped
   *  must not be implemented in engine class
   */
  explicit EntityPointer(const impl_t impl) : impl_(impl) {}

    
  /**
   * \brief	Generic assignment operator: Any entity pointer can be assigned
   * 			to this one.
   * \tparam EP_IMPL2 The underlying implementation type of the other 
   * 					EntityPointer wrapper.
   * \param	rhs	The entity pointer which should be assigned.
   * \return	A copy of "this" object.
   * \remarks The type EP_IMPL2 must inherit from  the
   * 			gridTraits::entityPointer_t<codim>() class!
   */
  template<class EP_IMPL2>
  EntityPointer<gridTraits_t,impl_t>& 
    operator=(const EntityPointer<gridTraits_t, EP_IMPL2>& rhs) {
      static_assert(std::is_reference<decltype(impl_.operator=(rhs.impl_))>::value, "operator=() in implementation must return a reference!");
      static_assert(std::is_base_of<typename gridTraits_t::template entityPointer_t<codim>,EP_IMPL2>::value,
        "template type EP_IMPL2 must inherit from gridTraits_t::entityPointer_t<codimension>");
      impl_.assign(rhs.impl_);
      return *this;
  }

  // Redirect from default assignment operator to template version
  // cf. (http://stackoverflow.com/questions/5625790/template-assignment-operator-overloading-mystery)
  EntityPointer<gridTraits_t,impl_t>& operator=(const EntityPointer<gridTraits_t,impl_t>& rhs) {
    return operator=<impl_t>(rhs);
  }


  /**
   * \brief	Indirection operator.
   * \remarks The underlying implementation must provide the method dereference()
   * 			which is called by this wrapper.
   */
  const Entity<gridTraits_t,codim>& operator*() const {
    static_assert(std::is_reference<decltype(impl_.dereference())>::value, "dereference() in implementation must return a reference!");
    return impl_.dereference();
  }

  /**
   * \brief	Member dereference operator.
   * \remarks The underlying implementation must provide the method dereference()
   * 			which is called by this wrapper.
   * \remarks We return only a const reference because EntityPointers are
   * 			immutable + if the EntityPointer is initialized by an Entity,
   * 			we only get a const reference to construct the pointer from...
   */
  const Entity<gridTraits_t,codim>* operator->() const {
    static_assert(std::is_reference<decltype(impl_.dereference())>::value, "dereference() in implementation must return a reference!");
    return &impl_.dereference();
  }

  /**
   * \brief	Do the two EntityPointers point to the same entity?
   * \param	rhs	The other entity pointer.
   * \remarks This is a generic operator so that we can compare e.g.
   * 			a LeafIterator to a LevelIterator.
   */
  template<class EP_IMPL2>
  bool operator==(const EntityPointer<gridTraits_t,EP_IMPL2>& rhs) const {
    static_assert(std::is_base_of<typename gridTraits_t::template entityPointer_t<codim>,EP_IMPL2>::value,
      "template type EP_IMPL2 must inherit from gridTraits_t::entityPointer_t<codimension>");
    return (impl_ == rhs.impl_);
  }

  /**
   * \brief	Generic Inequality operator.
   * \return	true if the two EntityPointers do not point to the same entity.
   * \remarks This function must not be implemented in the underlying engine
   * 			class. Naturally the negation of operator== is used.
   */
  template<class EP_IMPL2>
  bool operator!=(const EntityPointer<gridTraits_t,EP_IMPL2>& rhs) const {
    static_assert(std::is_base_of<typename gridTraits_t::template entityPointer_t<codim>,EP_IMPL2>::value,
                  "template type EP_IMPL2 must inherit from gridTraits_t::entityPointer_t<codimension>");
    return !(impl_ == rhs.impl_);
  }

  /// On which refinement level does the entity lie to which the pointer points?
  // This method is redundant and is only here for efficiency...
  size_type level() const {
    return impl_.level();
  }

protected:
  /// The underlying implementation...
  impl_t impl_;


private:
  
};

#endif // entityPointer_HPP__
