/********************************************************************
  created:	2013/06/19
  created:	19:6:2013   16:06
  filename: 	indexSet.hpp
  author:		Raffael Casagrande
  
  purpose:	Provide abstract interface for an IndexSet
*********************************************************************/


#ifndef indexSet_HPP__
#define indexSet_HPP__

#include <set>

namespace eth{
  namespace grid{

    /**
     * \brief	Assigns a consecutive, zero based index which is unique for
     * 			all boundary intersections and all entities of the same RefElType.
     * 			This index becomes invalid when the grid changes.
     * \tparam GRID_TRAITS	a class which models \ref GridTraitsDoc .
     * \tparam IS_IMPL		The underlying implementation class, must be
     * 						a model of IndexSet. 
     * 
     * 
     *
     * IndexSets provide a map,
     * \f[ m : E \mapsto \mathbb{N}^+ \f]
     * where \f$ E \f$ is a subset of all entities in the mesh (e.g. all 
     * Leaf Entities or all Entities on a given level) and 
     * \f$\mathbb{N}^+\f$ is the set of positive natural numbers (including 0).
     * 
     * \remark The subset \f$ E \f$ is defined in the implementation of the
     * 		   Grid. However \e every Grid provides a Leaf and Level Index
     * 		   Set (Grid::leafIndexSet() and Grid::levelIndexSet(), 
     * 		   respectively GridView::indexSet() )
     * 		   
     * Let us introduce the subset
     * \f[ E^t = \{ e \in E | \,  \text{ref El Type}(e)=t \}. \f]
     *  
     *
     * Index Sets have two unique properties which must be respected by the
     * underlying implementation:
     * - The index assigned to an entity is unique within \f$ E^t \f$  
     *	 , i.e. for all   
     *   \f$ e,e' \in E^t \f$ such that \f$ e \neq e' \f$ we have 
     *   \f$ m(e) \neq m(e') \f$.
     * - The index for all entities of the same codimension and reference
     *   element type is consecutive and zero based, i.e. for all   
     *   \f$ e \in E^t \f$ we have \f$ 0\leq m(e) < |E^t| \f$.
     *   	 
     * \note The mapping function of IndexSets has constant complexity but
     * 		 this comes at the price that the IndexSet becomes invalid when
     * 		 the mesh changes (e.g. is refined).
     * \note The indices to which the entities are mapped are arbitrary
     * 		 (dependent on the implementation).
     * 
     * Usually IndexSet 's are not used directly, but in conjunction with a
     * mapper. Mappers are automatically provided by the Interface and must
     * \e not be implemented by the Grid-implementor. Using mappers, the user
     * can map entities of different codimension/referenceElementTypes to a
     * unique, consecutive index (i.e. the mapped value of 
     * \f$ e_1 \in E^{t_1} \f$ and \f$ e_2 \in E^{t_2} \f$
     * is not the same for all \f$ t_1,t_2\f$ as long as 
     * \f$ e_1 \neq e_2 \f$).
     * 
     *
     *
     * <H3>Notes for Grid-Implementors </H3>
     * The IndexSet uses the CRTP pattern much in the same way as the Entity
     * class does. The reason for this is that IndexSet s cannot be copied
     * are assigned to. Therefore CRTP can be used and is more convenient 
     * for the implementor because he doesn't have to care about maintaining
     * an engine wrapper class.
     * 
     * 
     */
    template<class GRID_TRAITS, class IS_IMPL>
    class IndexSet {
    public:
      typedef GRID_TRAITS gridTraits_t;
      typedef IS_IMPL impl_t;
      typedef typename gridTraits_t::size_type size_type;

      /// Forbid assignment
      IndexSet& operator=(const IndexSet& rhs) = delete;


      /**
       * \brief	Map an entity of given codimension to the index.
       * \tparam	CODIM the codimension of the entity which should be
       * 				  mapped.
       * \param	e	The entity (codimension=CODIM) which should be mapped.
       * \warning	The behaviour of this method is undefined if the 
       * 			the domain \f$ E \f$ of the IndexSet does not contain e.
       * 			(i.e. if contains(e)==false)
       */
      template<int CODIM>
      size_type index(const Entity<gridTraits_t,CODIM>& e) const;


      /**
       * \brief Map a subEntity of a given element (codim=0)
       *  of a given codimension to it's index.
       * \tparam CODIM The codim of the subentity
       * \param e     The element from which a subentity should be mapped (codim=0).
       * \param i     The index of the subEntity.
       * \return      index of the specified subentity.
       *  
       *  \remark  This method produces exactly the same as 
       *    \code
       *           index(e.subEntity&lt;codim&gt;(i))
       *    \endcode (not perfectly legal code). 
       *    The current method is only here as a shortcut and to offer
       *    (possibly) superior performance.
       */
      template<int CODIM>
      size_type subIndex(const Entity<gridTraits_t,0>& e, 
         size_type i) const;

      /**
       * \brief	Gets a list of all Reference Elements 
       * 			which are present in the domain
       * 			\f$ E \f$ of the IndexSet.
       */
      std::set<eth::base::RefElType> refElTypes() const;

      /**
       * \brief	Get total number of entities of a given ReferenceElement
       * 			type in the entity set \f$ E \f$.
       * \param	type	The type reference element type.
       * \remark	If there are no Reference elements of the given types,
       * 			return 0.
       */
      size_type size(eth::base::RefElType type) const;

      /// Is the entity e contained in the set \f$ E \f$?
      template<int CODIM>
      bool contains(const Entity<gridTraits_t,CODIM>& e) const;

      /**
       * \brief	Is the given subEntity contained in \f$ E \f$ ?
       * \tparam	CODIM	Codimension of the subEntity.
       * \param	e	 	The element of which we take the subEntity
       * \param	i	 	Zero-based index of the subEntity
       * 
       */
      template<int CODIM>
      bool subContains(const Entity<gridTraits_t,0>& e, size_type i) const;

    protected:
      /// Forbid copy at least on interface level.
      IndexSet(const IndexSet& other);
      /// Default constructor
      IndexSet();
      /// Destructor
      ~IndexSet() {}

    private:
      /// CRTP trick
      impl_t& asImp() { return static_cast<impl_t&> (*this);}
      const impl_t& asImp()const {return static_cast<const impl_t&>(*this);}
    };
  }
}

// Implementation
//////////////////////////////////////////////////////////////////////////
#include "eth_base/CRTP_check.hpp"

namespace eth{
  namespace grid{

    template<class GRID_TRAITS, class IS_IMPL>
    IndexSet<GRID_TRAITS,IS_IMPL>::IndexSet() {}

    template<class GRID_TRAITS,class IS_IMPL>
    IndexSet<GRID_TRAITS,IS_IMPL>::IndexSet(const IndexSet& /* other */) {  }

    template<class GRID_TRAITS, class IS_IMPL>
    template<int CODIM>
    typename IndexSet<GRID_TRAITS,IS_IMPL>::size_type
      IndexSet<GRID_TRAITS,IS_IMPL>::index
      (const Entity<gridTraits_t,CODIM>& e) const {
        CHECK_CRTP_RETURN(asImp().index(e));
    }

    template<class GRID_TRAITS, class IS_IMPL>
    template<int CODIM>
    typename IndexSet<GRID_TRAITS,IS_IMPL>::size_type
      IndexSet<GRID_TRAITS,IS_IMPL>::subIndex(const Entity<gridTraits_t,0>& e, 
      size_type i) const {
        CHECK_CRTP_RETURN(asImp().template subIndex<CODIM>(e,i));
    }

    template<class GRID_TRAITS, class IS_IMPL>
    std::set<eth::base::RefElType> 
      IndexSet<GRID_TRAITS,IS_IMPL>::refElTypes() const {
      CHECK_CRTP_RETURN(asImp().refElTypes());
    }

    template<class GRID_TRAITS, class IS_IMPL>
    typename IndexSet<GRID_TRAITS, IS_IMPL>::size_type 
      IndexSet<GRID_TRAITS, IS_IMPL>::size(eth::base::RefElType type) const {
      CHECK_CRTP_RETURN(asImp().size(type));
    }

    template<class GRID_TRAITS, class IS_IMPL>
    template<int CODIM>
    bool IndexSet<GRID_TRAITS, IS_IMPL>::contains
      (const Entity<gridTraits_t,CODIM>& e) const {
      CHECK_CRTP_RETURN(asImp().contains(e));
    }

    template<class GRID_TRAITS, class IS_IMPL>
    template<int CODIM>
    bool IndexSet<GRID_TRAITS,IS_IMPL>::subContains
      (const Entity<gridTraits_t,0>& e, size_type i) const {
        CHECK_CRTP_RETURN(asImp().template subContains<CODIM>(e,i));
    }


  }
}




#endif // indexSet_HPP__
