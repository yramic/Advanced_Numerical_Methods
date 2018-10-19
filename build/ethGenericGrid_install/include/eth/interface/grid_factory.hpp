/********************************************************************
  created:	2013/07/26
  created:	26:7:2013   17:00
  filename: 	grid_factory.hpp
  author:		Raffael Casagrande
  
  purpose:	Provide an abstract interface for constructing unstructured
            grids...
*********************************************************************/
//          Copyright Raffael Casagrande 2013 - 2015.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#ifndef _HPP_2D7C59D5_868A_4972_B00D_33EBD84DAB90
#define _HPP_2D7C59D5_868A_4972_B00D_33EBD84DAB90


namespace eth {
  namespace grid{
    template<class GRID>
    class GridFactory;
  }
}

// Implementation
//////////////////////////////////////////////////////////////////////////

#include "eth_base.hpp"
#include <string>
#include <boost/lexical_cast.hpp>
#include <boost/function.hpp>

namespace eth {
  namespace grid {

    /**
     * \brief Provide an implementation independent way of creating the
     *        unstructured, macro grid.
     * 
     * Generally every Grid Implementation provides a grid factory which can
     * then in turn be used to initialize a grid. The present interface 
     * specifies the minimal functionality which should be supported by the
     * underlying grid class. 
     * 
     * The grid creation process can be split into five stages:
     * 1) You insert all the nodes which are necessary to define the geometry of  
     *    the elements. (You can also insert auxiliary nodes which are needed
     *    to define higher order elements.)
     * 2) Insert all the elements which are defined by their GeometryType and  
     *    indices of the nodes.
     * 3) (optional) insert BoundarySegments.
     * 4) Call CreateGrid() to finalize grid creation.  
     * 5) To determine the insertion index of grid nodes,elements and boundary  
     *    intersections you can now call insertionIndex().
     *    
     * In principal the first three steps can also occur in a mixed fashion,
     * the only requirement is that the nodes needed to define an element or
     * boundarySegment must be introduced before the element/boundarySegment.
     * 
     * In order to allow the user of the GridInterface to attach data to the
     * Grids entities the concept of an insertion index is needed:
     * The grid implementation guarantees that all nodes,elements and
     * boundary segments which are inserted through the GridFactory receive a
     * unique, zero-based and consecutive insertion index. I.e. if you insert
     * N nodes, they will be numbered according to the insertion order, the
     * first one will have insertion index 0. (A node can still have the same
     * insertion index as a Boundary segment or an element!).
     *  
     * Although this is redundant, the insertionIndices are also returned by the
     * insertNode(), insertElement() and insertBoundarySegment() methods.
     * 
     * After the grid has been created the user can ask for the insertion indices
     * of nodes, elements and boundary interfaces. This is useful to attach e.g.
     * material data to elements/nodes or to assign Dirichlet boundary conditions
     * to boundary Segments.
     *
     */
    template <class GRID_TRAITS>
    class GridFactory {
    public:
      typedef GRID_TRAITS gridTraits_t;
      /// Underlying implementation class.
      typedef typename gridTraits_t::gridFactory_t impl_t;
      typedef typename gridTraits_t::size_type size_type;
      /// Standard fixed size matrix type:
      template<int ROWS, int COLS>
      using matrix_t = typename gridTraits_t:: template fixedSizeMatrix_t<ROWS,COLS>;

      static const int dimWorld = gridTraits_t::dimWorld;
      static const int dimMesh = gridTraits_t::dimMesh;

      /// forbid copy and assignment (CRTP!)
      GridFactory(const GridFactory&) = delete;
      GridFactory& operator=(const GridFactory&) = delete;


      /// @{
      /// Functions which shall be called before createGrid()

      /**
       * \brief Insert a fixed number of nodes into the coarse grid.
       * \param coords  The coordinates of the node as a column vector.
       * \return The insertion index of the inserted point.
       * \note The nodes which are introduced here must not necessarily be
       *       "main nodes", i.e. nodes which belong to the corners of a
       *       reference element. They can also be "auxiliary nodes" which are
       *       used to define higher order elements. The GridFactory is
       *       responsible for distinguishing between Main and auxilliary nodes.
       * \note Every node which you insert using this routine gets a zero-based,
       *       consecutive index which can then be used with the
       *       insertElement() and insertBoundarySegment() routines.
       */
      size_type insertNode
        (const typename gridTraits_t::template fixedSizeMatrix_t<gridTraits_t::dimWorld,1>& coords) {
        CHECK_CRTP_RETURN(asImp().insertNode(coords));
      }


      /**
       * \brief Inserts an element (codim=0) which is defined by the given nodes
       *        (specified through their insertion index)
       * \param type      The GeometryType of the element which is inserted.
       * \param nodes     The insertion indices (from insertNodes()) of the nodes
       *                  which are needed to define the geometry of this element.
       * \return The insertion index of the inserted element.
       * \note Every element which is inserted with this routine gets automatically
       *       a zero-based consecutive index which can later on be used to 
       *       attach data to the grid.
       * \note Make sure the element is not inverted!
       */
      size_type insertElement(const eth::base::GeometryType& type, 
        const std::vector<size_type>& nodes) {
          ETH_ASSERT_MSG(nodes.size()==eth::base::GeometryTypeInfos::numNodes(type),
            "nodes.size() does not match the number of nodes of the geometry type.");
          ETH_ASSERT_MSG(eth::base::GeometryTypeInfos::dimension(type) == 
            gridTraits_t::dimMesh,
            (std::string("Geometry Type ") + eth::base::getGeometryTypeName(type) +
            std::string(" has dimension ") + 
            boost::lexical_cast<std::string>(eth::base::GeometryTypeInfos::dimension(type)) +
            std::string(" but the mesh has dimension ") +
            boost::lexical_cast<std::string>(gridTraits_t::dimMesh) +
            std::string(". This is therefore not an element.")).c_str());

         CHECK_CRTP_RETURN(asImp().insertElement(type,nodes));
      }


      /**
       * \brief Inserts a boundary segment which is defined by its nodes. This
       *        boundary segment should be a face of an element!
       * \param type  GeometryType of the face.
       * \param nodes The nodes which define the face.
       * \return  The insertion index of the boundarySegment.
       * \note  There is no need to call this function for all boundary segments
       *        of the grid. The Grid implementation will detect them 
       *        automatically. Nevertheless this function is useful because it
       *        assigns explicitly an insertion index to the boundary segment which can
       *        later on be used to get a reference to the face and attach e.g.
       *        Dirichlet data to it. The boundary faces which are not
       *        inserted explicitly and which are detected implicitly by the
       *        GridFactory receive an arbitrary insertion index.
       */
      size_type insertBoundarySegment(const eth::base::GeometryType type,
        const std::vector<size_type>& nodes) {
          ETH_ASSERT_MSG(nodes.size()==eth::base::GeometryTypeInfos::numNodes(type),
            "nodes.size() does not match the number of nodes of the geometry type.");
          ETH_ASSERT_MSG(eth::base::GeometryTypeInfos::dimension(type)==dimMesh-1,
            "The supplied geometry type cannot be a face of the mesh.");
          CHECK_CRTP_RETURN(asImp().insertBoundarySegment(type,nodes));
      }

      /**
       * \brief Inserts a boundary segment with arbitrary parameterization. 
       * \param type            The type of the face from which we map.
       * \param nodes           The nodes characterizing the face from which we
       *                        map.
       * \param parametrization A function object which takes a vector of length
       *                        dimMesh-1 and returns a vector of length dimWorld.
       *                        This mapping effectively describes the 
       *                        boundary parameterization.
       * \return The insertion index of this boundary segment.
       * \note Each boundary segment which is inserted by this method or by
       *       insertBoundarySegment(const eth::base::GeometryType type, const std::vector<size_type>& nodes)
       *       receives a unique insertion index which can be used after grid
       *       creation to identify a boundary segment.
       *       
       */
      size_type insertBoundarySegment(const eth::base::GeometryType type, 
        const std::vector<size_type>& nodes, const boost::function<
        matrix_t<dimWorld,1>(const matrix_t<dimMesh-1,1>&)>& parametrization) {
          ETH_ASSERT_MSG(nodes.size()==eth::base::GeometryTypeInfos::numNodes(type),
            "nodes.size() does not match the number of nodes of the geometry type.");
          ETH_ASSERT_MSG(eth::base::GeometryTypeInfos::dimension(type)==dimMesh-1,
            "The supplied geometry type is not a face of the mesh.");
          CHECK_CRTP(asImp().insertBoundarySegment(type,nodes,parametrization));
      }

      /// @}

      /**
       * \brief Creates an instance of the grid and passes it over as a pointer.
       *        Call this method once you have inserted all necessary nodes,
       *        elements and boundary segments.
       * \return  nullptr if it fails, else the new grid.
       * \warning The receiver takes the responsibility to deallocate the 
       *          grid.
       */
      typename gridTraits_t::gridImpl_t* createGrid() {
        CHECK_CRTP_RETURN(asImp().createGrid());
      }

      
      ///@{ Methods which can be called after createGrid() has been called.

      /**
       * \brief Retrieve the insertion index of an element.
       *        
       * This method can be used when the Grid has been create through
       * createGrid(). It is useful e.g. to associate material parameters with
       * elements.
       * 
       * \note It is also possible to pass elements which result from the refinement
       *       of other elements. In this case the insertionIndex of the 
       *       corresponding macro grid element is returned.
       * 
       */
      size_type insertionIndex(const Entity<gridTraits_t,0>& entity) const {
        CHECK_CRTP_RETURN(insertionIndex(entity));
      }

      /**
       * \brief Retrieve the insertion index for a node.
       * 
       * \note It is not possible to retrieve the index of an auxilliary node,
       *       because only main nodes are entities.
       */
      size_type insertionIndex(const Entity<gridTraits_t,dimMesh>& entity) const {
        CHECK_CRTP_RETURN(insertionIndex(entity));
      }

      /**
       * \brief Retrieve insertion index of a boundary intersection. 
       * 
       * \note If an intersection with boundary() == false is passed,
       *       an error is thrown.
       * \note It is also possible to pass intersections which result from refining
       *       the macro grids intersections. In this case the insertionIndex
       *       of the corresponding macro grid intersection is returned.
       */
      size_type insertionIndex(const Intersection<gridTraits_t>& intersection) const {
        ETH_ASSERT_MSG(intersection.boundary(),"This is not a boundary intersection.");
        CHECK_CRTP_RETURN(insertionIndex(intersection));
      }

      /// @}


    protected:
      /// Hide constructor (CRTP!)
      explicit GridFactory() {}
      /// This object cannot be destructed directly (CRTP) you must cast it to the actual implementation.
      ~GridFactory() {}

      

    private:
      /// CRTP trick
      impl_t& asImp() { return static_cast<impl_t&> (*this);}
      const impl_t& asImp()const {return static_cast<const impl_t&>(*this);}
    };
  }
}

#endif // _HPP_2D7C59D5_868A_4972_B00D_33EBD84DAB90