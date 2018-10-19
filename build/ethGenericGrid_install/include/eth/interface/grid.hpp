/********************************************************************
  created:	2013/06/17
  created:	17:6:2013   14:31
  filename: 	C:\Users\chracas\Documents\ETHGrid\Libs\interface\grid.hpp
  file path:	C:\Users\chracas\Documents\ETHGrid\Libs\interface
  author:		Raffael Casagrande
  
  purpose:	Provide an abstract interface for a generic grid class
*********************************************************************/

#ifndef grid_HPP__
#define grid_HPP__

// Declaration
namespace eth{
  namespace grid{
    template<class GRID_TRAITS>
    class Grid;
  }
}

#include "entity.hpp"
#include "grid_view.hpp"
#include "entity_iterator.hpp"
#include "entity_collection.hpp"
#include "entity_pointer.hpp"
#include "entity_iterator.hpp"
#include "entity_seed.hpp"
#include "id_set.hpp"
#include "index_set.hpp"
#include "grid_view_types.hpp"
#include "eth_base/CRTP_check.hpp"
#include "boundary_index_set.hpp"
#include <boost/function.hpp>

#include <type_traits>


/**
 * \brief	Specifies the abstract interface for any Grid class.
 * 
 *
 * \tparam	GridTraits	A Traits class which specifies all types appearing in
 * 						the actual implementation (see \ref GridTraitsDoc ).
 * \remarks The Grid class implements the CRTP pattern. This is needed because
 * 			usually methods in the implementation need access to the this pointer.
 */
template<class GRID_TRAITS>
class eth::grid::Grid {
public:
  /// the traits class which describes the generic Grid Implementation
  typedef GRID_TRAITS gridTraits_t;
  /// The underlying implementation (CRTP)
  typedef typename gridTraits_t::gridImpl_t impl_t;
  /// Type which is used to count something (see \ref GridTraitsDoc )
  typedef typename gridTraits_t::size_type size_type;

  // forbid copy and assignment:
  Grid(const Grid<gridTraits_t>&) = delete;
  Grid<gridTraits_t>& operator=(const Grid<gridTraits_t>&) = delete;

  /**
   * \brief	Gets the maximum refinement level of this grid
   * \remarks if refinement is not supported this returns 0.
   * \return	the maximum possible refinement level.
   */
  size_type maxLevel() const {
    CHECK_CRTP_RETURN(asImp().maxLevel())
  }

  
  /**
   * \brief	The number of entities on a given level of the given codimension.
   * \param	level	The level.
   * \param	codim	The codim.
   * \return	The number of entities on a given level of the given codim.
   */
  size_type size(size_type level, int codim) const {
    CHECK_CRTP_RETURN(asImp().size(level,codim));
  }

  
  /**
   * \brief	The number of leaf entities of given codimension.
   * \param	codim	The codim.
   * \return	The number of leaf entities of given codimension.
   */
  size_type size(int codim) const  {
    CHECK_CRTP_RETURN(asImp().size(codim));
  }

  /// How many entities of given Reference Element type do there exist in total?
  size_type size(eth::base::RefElType refElType) const {
    CHECK_CRTP_RETURN(asImp().size(refElType));
  }

  /// Number of entities with given ReferenceElement type on given level?
  size_type size(size_type level, eth::base::RefElType refElType) const {
    CHECK_CRTP_RETURN(asImp().size(level,refElType));
  }



  /**
   * \brief	Return a level view on the given level.
   * \param	level	The refinement level.
   * \return	A GridView object representing the level grid view.
   */
  typename eth::grid::GridView<typename gridTraits_t::template viewTraits_t<GridViewTypes::LevelView>> levelView(size_type level) const {
    CHECK_CRTP_RETURN(asImp().levelView(level));
  }
  
  /**
   * \brief	Return the leaf view of all entities.
   * \return	A GridView object representing all the leaf entities.
   */
  eth::grid::GridView<typename gridTraits_t::template viewTraits_t<GridViewTypes::LeafView>> leafView() const {
    CHECK_CRTP_RETURN(asImp().leafView());
  }


  ///**
  // * \brief	Returns all entities on a given level of given codimension.
  // * \tparam CODIM the codimension of the entities to return.
  // * \param	level	The refinement level
  // * \return	A collection of entities (through EntityCollection)
  // * \remarks This returns a handy entity collection which can be used with
  // * 			the range based for loops of C++11.
  // */
  template<int CODIM>
  const eth::grid::EntityCollection<EntityIterator<gridTraits_t,
    typename gridTraits_t::template viewTraits_t<GridViewTypes::LevelView>::template entityIterator_t<CODIM>>>
    levelEntities(size_type level) const {
      CHECK_CRTP_RETURN(asImp().template levelEntities<CODIM>(level));
  }

  /**
   * \brief	Returns all leaf entities of a given codimension.
   * \tparam CODIM the codimension of the entities to return
   * \return	A collection of entities (through EntityCollection)
   * \remarks This returns a handy entity collection which can be used with
   * 			the range based for loops of C++11.
   */
  template<int CODIM>
  const eth::grid::EntityCollection<EntityIterator<gridTraits_t,
        typename gridTraits_t::template viewTraits_t<GridViewTypes::LeafView>::template entityIterator_t<CODIM>>>
    leafEntities() const {
    CHECK_CRTP_RETURN(asImp().template leafEntities<CODIM>());
  }

  /**
   * \brief access the refinement class
   * \return a reference to the refinement class which specifies the refinement algorithm 
   */
  const typename gridTraits_t::refinementMethod_t& refinement() const {
    CHECK_CRTP_RETURN(asImp().refinement());
  }

  /**
   * \brief access the refinement class
   * \return a reference to the refinement class which specifies the refinement algorithm 
   */
  typename gridTraits_t::refinementMethod_t& refinement() {
    CHECK_CRTP_RETURN(asImp().refinement());
  }

  /**
   * \brief	Refine the grid refCount times using default refinement rule.
   * \param	refCount How many times should the grid be refined?
   */
  void globalRefine(size_type refCount) {
    CHECK_CRTP(asImp().globalRefine(refCount));
  }

  /**
   * \brief	Marks an entity (codim=0) to be refined/coarsened in a
   * 			subsequent adapt.
   * \param	refCount	Number of subdivisions that should be applied.
   * 						Negative value means coarsening.
   * \param	e			Entity that should be marked.
   * \return	true if it succeeds, false if it fails.
   */
  bool mark(size_type refCount, const Entity<gridTraits_t,0>& e ) {
    CHECK_CRTP_RETURN(asImp().mark(refCount, e));
  }

  /**
   * \brief	Returns how many times the given element (codim=0) is to be
   * 			subdivided in the subsequent adapt. (Zero if it is not refined)
   * \param	e	The entity for which the number of subdivisions should be 
   * 				retrieved
   * \return	number of times Entity e is subdivided.
   */
  size_type getMark(const Entity<gridTraits_t,0>& e) const {
    CHECK_CRTP_RETURN(asImp().getMark(e));
  }

  /**
   * \brief	Must be called after entities have been marked but before
   * 			adapt() is called. --> Sets the mightVanish flags of the
   * 			elements for the next adapt call.
   * \return	true if an entity may be coarsened during a subsequent adapt(),
   * 			false otherwise.
   */
  bool preAdapt() {
    CHECK_CRTP_RETURN(asImp().preAdapt());
  }

  
  /**
   * \brief	Refines respectively coarsens all marked entities if possible.
   * \return	true if at least one entity was refined.
   * 			
   * In order to refine a grid follow the following steps:
   * -# Mark elements (codim=0) with the mark() method.
   * -# call preAdapt() to set mightVanish flags on entity
   *		- if preAdapt() returns true you might have to save the solution  
   * -# call adapt() 
   *		- if adapt() return true, you might have to interpolate the (saved)  
   *			solution
   * -# call postAdapt() to remove the isNew flags from the elements.
   */
  bool adapt() {
    CHECK_CRTP_RETURN(asImp().adapt());
  }

  /**
   * \brief	Call this after adapt() has been called. (removes the isNew
   * 			flags from elements)
   */
  void postAdapt() {
    CHECK_CRTP(asImp().postAdapt());
  }

  
  /**
   * \brief	Get EntityPointer from entitySeed.
   * \param	seed	The entity seed from which we want to restore the entity.
   * \return	Entity pointer which points to the entity from which the seed
   * 			was created
   * \remarks An Entity seed is a low-memory footprint reference to an Entity.
   * 			It can be created by calling Entity::seed() and can be changed
   * 			back into an Entity by calling this method.
   */
  template<int CODIM>
  EntityPointer<gridTraits_t,
                typename gridTraits_t::template entityPointer_t<CODIM>>  
    entityPointer(const EntitySeed<gridTraits_t, CODIM>& seed ) const {
    CHECK_CRTP_RETURN(asImp().entityPointer(seed));
  }


  /**
   * \brief	return const reference to the grids id set.
   * \return	An Id Set.
   * \remarks Use this method to obtain an idSet which in turn maps an entity
   * 			uniquely to an index/id. (Mainly used by mappers)
   */
  const IdSet<gridTraits_t>& idSet() const {
    static_assert(std::is_reference<decltype(asImp().idSet())>::value, "idSet() in implementation must return a reference!");
    CHECK_CRTP_RETURN(asImp().idSet());
  }

  /**
   * \brief	obtain const reference to a Index set mapping all entities on a
   * 			given level to a unique index. (i.e. the IndexSet only maps
   * 			from a Subset of all entities to an index)
   * \param	level	The refinement level.
   * \return	An IndexSet
   */
  const IndexSet<gridTraits_t,
                 typename gridTraits_t::template viewTraits_t<GridViewTypes::LevelView>::indexSet_t> &
    levelIndexSet(size_type level) const {
    static_assert(std::is_reference<decltype(asImp().levelIndexSet(level))>::value, "levelIndexSet() in implementation must return a reference!");
    CHECK_CRTP_RETURN(asImp().levelIndexSet(level));
  }

  /**
   * \brief	obtain a const reference to a Index set mapping all leaf
   * 			entities to a unique index.
   * \return	An Index Set.
   */
  const IndexSet<gridTraits_t,
                 typename gridTraits_t::template viewTraits_t<GridViewTypes::LeafView>::indexSet_t> &
    leafIndexSet() const {
    static_assert(std::is_reference<decltype(asImp().leafIndexSet())>::value, "leafIndexSet() in implementation must return a reference!");
    CHECK_CRTP_RETURN(asImp().leafIndexSet());
  }

  /**
   * \brief obtain a const reference a an BoundaryIndexSet which maps all 
   *        boundary intersections (any level) to a unique index.
   */
  const BoundaryIndexSet<gridTraits_t> & boundaryIndexSet() const {
    static_assert(std::is_reference<decltype(asImp().boundaryIndexSet())>::value, 
      "boundaryIndexSet() in implementation must return a reference!");
    CHECK_CRTP_RETURN(asImp().boundaryIndexSet());
  }

  /**
   * \brief Register a callback which is called \e before the underlying
   *        grid changes (e.g. before refinement.)
   * \param listener  A listener object which will be called.
   * \return  A handle on this listener object which can later on be used
   *          to deregister the listener again (through 
   *          removePreUpdateListener())
   *          
   * The listener will be called before the grid changes, but when entities which
   * might vanish are known. (e.g. eth::grid::Entity< GRID_TRAITS, 0 >::mightVanish() has a meaning)
   *          
   * \warning If you the listener is bound to some member function of an object,
   *          it is absolutely necessary that listener is removed before
   *          the object is destroyed!
   */
  typename gridTraits_t::listenerHandle_t addPreUpdateListener(const boost::function<void()>& listener) {
    CHECK_CRTP_RETURN(asImp().addPreUpdateListener(listener));
  }

  /**
   * \brief Removes the pre update listener described by listenerHandle.
   * \param listenerHandle  Handle of the listener, is obtained by when
   *                        addPreUpdateListener() is called.
   */
  void removePreUpdateListener(const typename gridTraits_t::listenerHandle_t& listenerHandle) {
    CHECK_CRTP_RETURN(asImp().removePreUpdateListener(listenerHandle));
  }

  /**
   * \brief Register a callback which is called \e after the underlying
   *        grid changes (e.g. after refinement.)
   * \param listener  A listener object which will be called.
   * \return  A handle on this listener object which can later on be used
   *          to deregister the listener again (through 
   *          removePreUpdateListener())
   *          
   * The listener will be called after the grid has changed and the new entities
   * are known. (e.g. eth::grid::Entity< GRID_TRAITS, 0 >::isNew() has a meaning!)
   *          
   * \warning If you the listener is bound to some member function of an object,
   *          it is absolutely necessary that listener is removed before
   *          the object is destroyed!
   */
  typename gridTraits_t::listenerHandle_t addPostUpdateListener(const boost::function<void()>& listener) {
    CHECK_CRTP_RETURN(asImp().addPostUpdateListener(listener));
  }

  /**
   * \brief Removes the pre update listener described by listenerHandle and
   *        which was previously registered through addPostUpdateListener().
   * \param listenerHandle  Handle of the listener, is obtained by when
   *                        addPostUpdateListener() is called.
   */
  void removePostUpdateListener(const typename gridTraits_t::listenerHandle_t& listenerHandle) {
    CHECK_CRTP_RETURN(asImp().removePostUpdateListener(listenerHandle));
  }
  


protected:
  /// Hide constructor (Barton-Nackmann trick)
  explicit Grid() {}

  /// This object cannot be destructed directly (CRTP) you must cast it to the actual implementation.
  ~Grid() {}

private:
  /// CRTP trick
  impl_t& asImp() { return static_cast<impl_t&> (*this);}
  const impl_t& asImp()const {return static_cast<const impl_t&>(*this);}
};






#endif // grid_HPP__
