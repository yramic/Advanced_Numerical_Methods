/********************************************************************
  created:	2013/06/17
  created:	17:6:2013   17:11
  filename: 	C:\Users\chracas\Documents\ETHGrid\Libs\interface\gridView.hpp
  file path:	C:\Users\chracas\Documents\ETHGrid\Libs\interface
  author:		Raffael Casagrande
  
  purpose:	Provide an abstract interface for a GridView.
*********************************************************************/


#ifndef gridView_HPP__
#define gridView_HPP__

#include "grid_view_types.hpp"

namespace eth {
  namespace grid {
    template<class VIEW_TRAITS>
    class GridView;
  }
}

// Implementation
//////////////////////////////////////////////////////////////////////////

#include "index_set.hpp"
#include "intersection_iterator.hpp"

/**
 * \brief	Provides an abstract interface for a GridView. 
 * 			GridView's simplify the 
 * 			hierarchical Grid structure by flattening the hierarchy and
 * 			thus providing a simpler Interface which is often enough.
 * \tparam VIEW_TRAITS a class which models ViewTraits (see \ref ViewTraitsDoc )
 * 					   and determines the underlying implementation class.
 * 
 * There are two general types of GridViews:
 * - LeafGridViews represent all leaf Entities of the Grid  
 * - LevelGridViews which represent all Entities on a given refinement level.
 * 
 * This class provides the following functionality:
 * - Reference to an indexSet which maps from the subset of all entities in this GridView.  
 * - Access to all entities through entities()  (can be used with range-based for loops!)
 * - Access to all intersections of a given element (codim=0)  through intersections().
 * - Reference to the underlying grid through grid().  
 * 
 * \remarks It is assumed that the underlying implementation is rather lightweight
 * 			and that it thus provides a copy constructor. Whenever the GridView
 * 			is copied, the implementation is copied as well. (This also means that the liftime of the object passed in the constructor
 * 			can end before the GridView is destructed)
 * 			If the Grid Implementator doesn't 
 * \remarks The GridView provides read-only access to all entities and becomes
 * 			invalid after refinement.
 */
template<class VIEW_TRAITS>
class eth::grid::GridView {
public:
  /// View Traits (see \ref ViewTraitsDoc)
  typedef VIEW_TRAITS viewTraits_t;
  /// the traits class which describes the generic Grid Implementation
  typedef typename viewTraits_t::gridTraits_t gridTraits_t;
  /// The underlying implementation (engine pattern)
  typedef typename viewTraits_t::viewImpl_t impl_t;

  typedef typename gridTraits_t::size_type size_type;

  /// Implementation Constructor, provide implementation as argument
  explicit GridView(impl_t impl) : impl_(std::move(impl)) {

  }

  ///Copy constructor
  GridView(const GridView<viewTraits_t>& other)
    : impl_(other.impl_){

  }

  /**  
  *  \brief Assignment operator
  *  \remark Assignment is needed because the user must be able to update
  *  		   an already existing GridView object by a new one after refinement.
  *  		   (See \ref LifetimeDoc for more infos)
  */
  GridView<viewTraits_t>&  operator=(const GridView<viewTraits_t>& rhs) {
    static_assert(std::is_reference<decltype(impl_.operator=(rhs))>::value, "operator=() in implementation must return a reference!");
    impl_ = rhs.impl_;
    return *this;
  }

  /**
   * \brief	Return reference to the full grid structure
   */
  const Grid<gridTraits_t>& grid() const {
    static_assert(std::is_reference<decltype(impl_.grid())>::value, "grid() in implementation must return a reference!");
    return impl_.grid();
  }

  /**
   * \brief	Returns an index set which maps from the subset of all entities
   * 			represented by this GridView. (Entities on given level respectively
   * 			entities on leaf level)
   * \return	An IndexSet
   */
  const IndexSet<gridTraits_t,typename viewTraits_t::indexSet_t>& indexSet() const {
    static_assert(std::is_reference<decltype(impl_.indexSet())>::value, "indexSet() in implementation must return a reference!");
    return impl_.indexSet();
  }

  /**
   * \brief	return number of entities of given codimension
   * \param	codim	The codimension.
   * \return	Number of entities of given codim.
   */
   size_type size(int codim) const {
    return impl_.size(codim);
  }

   /// Number of entities in this GridView of the given ReferenceElement type.
   size_type size(eth::base::RefElType refElType) const {
     return impl_.size(refElType);
   }

  /**
   * \brief	Is the given entity contained in this GridView?
   * \param	e	The entity.
   * \return	true if e is contained in this GridView, otherwise false.
   */
  template<int CODIM>
  bool contains(const Entity<gridTraits_t, CODIM>& e) const {
    return impl_.contains(e);
  }

  /**
  * \brief	Returns a collection of all entities of given Codimension in this
  * 			grid view. (can be used with range-based for loops)
  * \tparam CODIM the codim of the entities through which we want to traverse.
  */
  template<int CODIM>
  const 
  EntityCollection<EntityIterator<typename viewTraits_t::gridTraits_t,
                                  typename viewTraits_t::template entityIterator_t<CODIM>>> entities() const
  {
      return impl_.template entities<CODIM>();
  }

  /**
   * \brief	Returns a collection of all Intersections which entity e
   * 			has with other entities and the boundary
   * 			(can be used with range-base for loops).
   * \param	e	The Entity for which the intersections are to be obtained
   * \remark if this is a Leaf GridView and if e is not a leaf, 
   * 		   intersections.begin() == intersections.end()
   */
  const EntityCollection<IntersectionIterator<viewTraits_t>> 
    intersections(const Entity<gridTraits_t,0>& e) const {
    return impl_.intersections(e);
  }


private:
  impl_t impl_;
};


#endif // gridView_HPP__
