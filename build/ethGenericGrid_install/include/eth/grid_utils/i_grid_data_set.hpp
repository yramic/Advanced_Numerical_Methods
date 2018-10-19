
#ifndef ETH_I_GRID_DATA_SET_HPP
#define ETH_I_GRID_DATA_SET_HPP

#include <array>
#include <memory>

#include "eth_base.hpp"
#include "interface.hpp"

namespace eth {
  namespace grids {
    namespace utils {

      /**
       * \brief Engine interface class for a GridDataSet which is the basic component
       *        to attach any data to a grid's entities.
       * \tparam GDS_TRAITS A traits class which defines the underlying implementation class.
       * \note This is only the interface class which is not instantiable. If you
       *       want to create a new GridDataSet, take a look at the GridDataSet 
       *       class.
       *  
       * A GridDataSet allows the user to attach any kind of data to a set of entities.
       * This set can contain entities of different codimensions and different ReferenceElementTypes.
       * Typically the user can store a number of data entries with each entity, this
       * number (often also called Multiplicity) is defined on a per RefElType 
       * basis and is also specified by the user.
       *       
       * \note A IGridDataSet object is neither copyable nor moveable. If you want to
       *       create a deep copy of GridDataSet you must cast it into the actual implementation
       *       (CRTP pattern) and call it's copy constructor.
       * 
       * <H3>Note to implementor </H3>
       * A GDS_TRAITS class must provide the following public types:
       * type                      | meaning
       * --------------------------|----------------------------------------------------
       * impl_t                    | The actual class which implements the IGridDataSet interface (CRTP pattern)
       * dataType_t                | The type of data which is attached to the entities.
       * viewTraits_t              | A class which models \ref ViewTraitsDoc and specifies the underlying grid data structure.
       * 
       * In addition a a method with the following syntax must exist:
       * \code
       * static constexpr size_type multiplicity(RefElType).
       * \endcode
       * It specifies the number of data entries per reference element type and is
       * used in the std::array<dataType_t,impl_t::template multiplicity<TYPE> >& eth::io::IGridDataSet< GDS_IMPL >::data	(	const eth::grid::Entity< gridTraits_t, gridTraits_t::dimMesh-eth::base::ReferenceElement< TYPE >::dimension > & 	e)		
       * method.
       */
      template<class GDS_TRAITS>
      class IGridDataSet {
      public:
        typedef GDS_TRAITS traits_t;
        /// underlying CRTP implementation
        typedef typename traits_t::impl_t impl_t;
        /// The type of data which is attached to the entities.
        typedef typename traits_t::dataType_t dataType_t;
        /// This list specifies all refElTypes to which this gridDataSet attaches data.
        typedef typename traits_t::refElList_t refElList_t;
        typedef typename traits_t::viewTraits_t viewTraits_t;
        typedef typename viewTraits_t::gridTraits_t gridTraits_t;
        typedef eth::base::unsigned_t size_type;
      

        /**
         * \brief Gets the multiplicity for a given ReferenceElement type
         *        (i.e. the number of data entries stored for this RefElType)
         * \param type  The reference element type.
         * \return  The multiplicity, resp. number of data entries stored with
         *          ReferenceElement type.
         * \note  The RefElType 's which are mapped by this GridDataSet are
         *        publicly accessible through the public refElList_t type.
         */
        static constexpr eth::base::unsigned_t getMultiplicity(eth::base::RefElType type) {
          return traits_t::multiplicity2(type);
        }


        /**
         * \brief read/write access to i-th data entry which is attached to entity e.
         * \param e             The const eth::grid::Entity&lt;gridTraits_t,CODIM&gt;
         *                      to which the data is attached.
         * \param i             return data stored in i-th data entry, 0<=i<multiplicity
         *                      where the multiplicity is defined through e.RefElType().
         * \return  the i-th data entry attached to entity e.
         */
        template<int CODIM>
        dataType_t& data(const eth::grid::Entity<gridTraits_t,CODIM>& e, size_type i) {
          CHECK_CRTP_RETURN(asImp().data(e, i));
        }

        /**
         * \brief read/write access to i-th data entry which is attached to sub-entity subIndex of entity e.
         * \tparam CODIM        The co-dimension of the sub-entity to which the data is attached
         * \param e             The const eth::grid::Entity&lt;gridTraits_t,0&gt;
         *                      to whose sub-entity the data is attached.
         * \param subIndex      The local index in the entity e of the sub-entity to which the data is attached.
         * \param i             return data stored in i-th data entry, 0<=i<multiplicity
         *                      where the multiplicity is defined through e.subEntity<CODIM>(subIndex).RefElType().
         * \return  the i-th data entry attached to entity e.
         *          
         *  
         *  \remark  This method produces exactly the same as 
         *    \code
         *           data(e.subEntity&lt;CODIM&gt;(subIndex),i)
         *    \endcode. 
         *    The current method is only here as a shortcut and to offer
         *    (possibly) superior performance.
         */
        template<int CODIM>
        dataType_t& subData(const eth::grid::Entity<gridTraits_t,0>& e, size_type subIndex, size_type i) {
          CHECK_CRTP_RETURN(asImp().template subData<CODIM>(e,subIndex,i));
        }

        /**
         * \brief read-only access to i-th data entry which is attached to entity e.
         * \param e             The const eth::grid::Entity&lt;gridTraits_t,CODIM&gt;
         *                      to which the data is attached.
         * \param i             return data stored in i-th data entry, 0<=i<multiplicity
         *                      where the multiplicity is defined through e.RefElType().
         * \return  the i-th data entry attached to entity e.
         */
        template<int CODIM>
        const dataType_t& data(const eth::grid::Entity<gridTraits_t,CODIM>& e, size_type i) const {
          CHECK_CRTP_RETURN(asImp().data(e, i));
        }

        /**
         * \brief read-only access to i-th data entry which is attached to sub-entity subIndex of entity e.
         * \tparam CODIM        The co-dimension of the sub-entity to which the data is attached
         * \param e             The const eth::grid::Entity&lt;gridTraits_t,0&gt;
         *                      to whose sub-entity the data is attached.
         * \param subIndex      The local index in the entity e of the sub-entity to which the data is attached.
         * \param i             return data stored in i-th data entry, 0<=i<multiplicity
         *                      where the multiplicity is defined through e.subEntity<CODIM>(subIndex).RefElType().
         * \return  the i-th data entry attached to entity e.
         *          
         *  
         *  \remark  This method produces exactly the same as 
         *    \code
         *           data(e.subEntity&lt;CODIM&gt;(subIndex),i)
         *    \endcode. 
         *    The current method is only here as a shortcut and to offer
         *    (possibly) superior performance.
         */
        template<int CODIM>
        const dataType_t& subData(const eth::grid::Entity<gridTraits_t,0>& e, size_type subIndex, size_type i) const {
          CHECK_CRTP_RETURN(asImp().template subData<CODIM>(e,subIndex,i));
        }

        /**
         * \brief read/write access to all data entries attached to entity e.
         * \tparam TYPE The ReferenceElement type of the entity which is being passed.
         * \param e The const eth::grid::Entity&lt;gridTraits_t,CODIM&gt;
         *          to which the data is attached.
         * \return  all data entries attached to entity e.
         * \note    The reference element type of entity `e` must be `TYPE`, i.e.
         *          `e.refElType() == TYPE`
         */
        template<eth::base::RefElType TYPE>
        std::array<dataType_t,traits_t::template multiplicity<TYPE>::val>& data
          (const eth::grid::Entity<gridTraits_t,gridTraits_t::dimMesh-eth::base::ReferenceElement<TYPE>::dimension>& e) {
            CHECK_CRTP_RETURN(asImp().template data<TYPE>(e));
        }

        /**
         * \brief read/write access to all data entries attached to sub-entity subIndex of entity e.
         * \tparam SUB_TYPE The ReferenceElement type of the sub-entity which is being requested.
         * \param e             The const eth::grid::Entity&lt;gridTraits_t,0&gt;
         *                      to whose sub-entity the data is attached.
         * \param subIndex      The local index in the entity e of the sub-entity to which the data is attached.
         * \return  all data entries attached to entity subentity e.
         *          
         *  \remark  This method produces exactly the same as 
         *    \code
         *           data<SUB_TYPE>(e.subEntity&lt;CODIM&gt;(subIndex))
         *    \endcode. (here, CODIM denotes the co-dimension of SUB_TYPE with respect to e)
         *    The current method is only here as a shortcut and to offer
         *    (possibly) superior performance.
         *    
         * \note    The reference element type of sub-entity `subIndex` of entity `e` must be `SUB_TYPE`, i.e.
         *          `e.subEntity<CODIM>(subIndex).refElType() == SUB_TYPE`.
         *          `CODIM` denotes the co-dimension of `SUB_TYPE` with respect to `e`.
         */
        template<eth::base::RefElType SUB_TYPE>
        std::array<dataType_t,traits_t::template multiplicity<SUB_TYPE>::val>& subData
          (const eth::grid::Entity<gridTraits_t,0>& e, size_type subIndex) {
            CHECK_CRTP_RETURN(asImp().template subData<SUB_TYPE>(e,subIndex));
        }

        /**
         * \brief read-only access to all data entries attached to entity e.
         * \tparam TYPE The ReferenceElement type of the entity which is being passed.
         * \param e The const eth::grid::Entity&lt;gridTraits_t,CODIM&gt;
         *          to which the data is attached.
         * \return  all data entries attached to entity e.
         * \note    The reference element type of entity `e` must be `TYPE`, i.e.
         *          `e.refElType() == TYPE`
         */
        template<eth::base::RefElType TYPE>
        const std::array<dataType_t,traits_t::template multiplicity<TYPE>::val>& data
          (const eth::grid::Entity<gridTraits_t,gridTraits_t::dimMesh-eth::base::ReferenceElement<TYPE>::dimension>& e) const {
            CHECK_CRTP_RETURN(asImp().template data<TYPE>(e));
        }

        /**
         * \brief read-only access to all data entries attached to sub-entity subIndex of entity e.
         * \tparam SUB_TYPE The ReferenceElement type of the sub-entity which is being requested.
         * \param e             The const eth::grid::Entity&lt;gridTraits_t,0&gt;
         *                      to whose sub-entity the data is attached.
         * \param subIndex      The local index in the entity e of the sub-entity to which the data is attached.
         * \return  all data entries attached to entity subentity e.
         *          
         *  \remark  This method produces exactly the same as 
         *    \code
         *           data<SUB_TYPE>(e.subEntity&lt;CODIM&gt;(subIndex))
         *    \endcode. (here, CODIM denotes the co-dimension of SUB_TYPE with respect to e)
         *    The current method is only here as a shortcut and to offer
         *    (possibly) superior performance.
         *    
         * \note    The reference element type of sub-entity `subIndex` of entity `e` must be `SUB_TYPE`, i.e.
         *          `e.subEntity<CODIM>(subIndex).refElType() == SUB_TYPE`.
         *          `CODIM` denotes the co-dimension of `SUB_TYPE` with respect to `e`.
         */
        template<eth::base::RefElType SUB_TYPE>
        const std::array<dataType_t,traits_t::template multiplicity<SUB_TYPE>::val>& subData
          (const eth::grid::Entity<gridTraits_t,0>& e, size_type subIndex) const {
            CHECK_CRTP_RETURN(asImp().template subData<SUB_TYPE>(e,subIndex));
        }

      protected:
        /// hide default constructor (CRTP)
        IGridDataSet() {}
        /// hide copy constructor (CRTP)
        IGridDataSet(const IGridDataSet& ) {}
        // hide destructor (CRTP)
        ~IGridDataSet() {}


      private:
        /// CRTP trick
        impl_t& asImp() { return static_cast<impl_t&> (*this);}
        /// CRTP trick
        const impl_t& asImp()const {return static_cast<const impl_t&>(*this);}
      };
    }
  }
}

#endif // ETH_I_GRID_DATA_SET_HPP

