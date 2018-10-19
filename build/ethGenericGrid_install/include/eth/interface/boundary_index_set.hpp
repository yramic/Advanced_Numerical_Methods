/********************************************************************
  created:	2013/09/03
  created:	3:9:2013   11:01
  filename: 	boundary_index_set.hpp
  author:		Raffael Casagrande
  
  purpose:	Abstract interface to describe a mapping from boundary
            faces to zero-based consecutive indices.
*********************************************************************/
//          Copyright Raffael Casagrande 2013 - 2015.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)


#ifndef _HPP_54FAEC64_F446_4B04_9DD7_446D4509F160
#define _HPP_54FAEC64_F446_4B04_9DD7_446D4509F160


// Declaration
namespace eth {
  namespace grid {
    template<class GRID_TRAITS>
    class BoundaryIndexSet;
  }
}

// Implementation
//////////////////////////////////////////////////////////////////////////
#include "intersection.hpp"
namespace eth {
  namespace grid {




    /**
    * \brief Abstract interface to describe a mapping from boundary
    *        faces to zero-based consecutive indices. This index 
    *        becomes invalid when the grid changes.
    * \tparam GRID_TRAITS a class which models \ref GridTraitsDoc . (this
    *                     determines the underlying implementation class)
    * 
    * A BoundaryIndexSet is very similar to an ordinary index set. The big 
    * difference is that BoundaryIndexSets map boundary intersections to an 
    * index. As in the case of IndexSets this index is also zero-based and
    * consecutive and it becomes invalid when the grid changes !
    * 
    * In comparison to IndexSets, BoundaryIndexSets work with any type of
    * boundary Intersection object. It does not matter to what level this 
    * intersection object belongs to. This makes sense because of the following
    * property:
    * 
    * The index of an arbitrary boundary intersection is always the same
    *            as the index of the corresponding boundary intersection in the 
    *            macro grid.*
    *            
    * Thus the maximum index is given by the number of boundary 
    * intersections in the macro grid.
    * 
    * Because BoundaryIndexSet s become invalid when the grid changes one
    * must combine it with an IdSet in a similar way as with the IndexSet.
    * 
    * \note The boundary index has nothing to do with the insertion index of
    *       a boundary segment which appears in class GridFactory.
    * 
    * <H3>Notes for Grid Implementors</H3>
    * The BoundaryIndexSet class uses the CRTP pattern in the same way as the
    * IndexSet class. This is more convenient for the implementor because
    * he doesn't have to care about maintaining an engine wrapper class.
    * 
    * \remark The BoundaryIndexSet functionality was not integrated into the
    *         IndexSet class because IndexSets are typically different for 
    *         different levels which is not true for boundary intersections.
    *         It was also not integrated into the Grid interface to facilitate
    *         the development of library to attach data to a grid.
    * 
    */
    template<class GRID_TRAITS>
    class BoundaryIndexSet {
    public:
      typedef GRID_TRAITS gridTraits_t;
      typedef typename gridTraits_t::boundaryIndexSet_t impl_t;
      typedef typename gridTraits_t::size_type size_type;

      // forbid assignment (CRTP pattern)
      BoundaryIndexSet& operator=(const BoundaryIndexSet& rhs) = delete;

      /**
       * \brief Returns a zero-based, consecutive index for a given boundary
       *        intersection in constant time.
       * \param i A boundary intersection (i.boundary==true).
       * \return  A zero-based, consecutive index which corresponds to the
       *          given intersection. The maximal index is given by
       *          Grid::numBoundarySegments().
       * 
       *  \remark see the documentation of the BoundaryIndexSet class for more details.
       */
      size_type index(const Intersection<gridTraits_t>& i) const {
        ETH_ASSERT_MSG(i.boundary(),"Can only map boundary intersections.");
        CHECK_CRTP_RETURN(asImp().index(i));
      }

      /**
       * \brief Returns the number of boundary intersections. size()-1 is
       *        therefore also the maximum value which is returned by the
       *        index() method.
       */
      size_type size() const {
        CHECK_CRTP_RETURN(asImp().size());
      }
      

    protected:
      // Forbid copy at least on interface level:
      BoundaryIndexSet(const BoundaryIndexSet& other);
      /// Default Constructor:
      BoundaryIndexSet() {}
      /// Destructor:
      ~BoundaryIndexSet(){}

    private:
      // CRTP trick...
      impl_t& asImp() { return static_cast<impl_t&> (*this);}
      const impl_t& asImp()const {return static_cast<const impl_t&>(*this);}
    };
  }
}





#endif // _HPP_54FAEC64_F446_4B04_9DD7_446D4509F160