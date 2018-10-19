
#ifndef ETH_GRID_DATA_SET_HPP
#define ETH_GRID_DATA_SET_HPP

#include <stdexcept>

#include <boost/lexical_cast.hpp>

#include "eth_base.hpp"
#include "interface.hpp"
#include "grid_view_factory.hpp"

namespace eth {
  namespace grids {
    namespace utils {

      // Forward declaration:
      template<class DATA_TYPE, class GRID_VIEW_FACTORY, class... MULTIPLICITIES>
      class GridDataSet;

      /** \brief struct to define at compile time the number of data entries to
       *         be stored in a grid data set for a given reference element
       *  \tparam TYPE the reference element for which data shall be stored
       *  \tparam MULTIPLICITY the number of data entries to be stored per
       *                       reference element of type TYPE
       */
      template<eth::base::RefElType TYPE, eth::base::unsigned_t MULTIPLICITY>
      struct MultiplicityPair{
        static constexpr eth::base::RefElType refElType = TYPE;
        static constexpr eth::base::unsigned_t multiplicity = MULTIPLICITY;
      };

      namespace detail {
        /// traits class for the GridDataSet class below.
        template<class DATA_TYPE, class GRID_VIEW_FACTORY, class... MULTIPLICITIES>
        struct GridDataSetTraits {
          static_assert(!std::is_same<DATA_TYPE,DATA_TYPE>::value,
            "The template parameter GRID_VIEW_FACTORY of class GridDataSet must be of type grids::utils::GridViewFactory.");
        };
        
        template<class DATA_TYPE, class GRID, grid::GridViewTypes GV_TYPE, class... MULTIPLICITIES> 
        struct GridDataSetTraits<DATA_TYPE, grids::utils::GridViewFactory<GRID,GV_TYPE>, MULTIPLICITIES...>
        {
          typedef DATA_TYPE dataType_t;
          typedef typename grids::utils::GridViewFactory<GRID,GV_TYPE>::viewTraits_t viewTraits_t;
          typedef GridDataSet<dataType_t,GridViewFactory<GRID,GV_TYPE>,MULTIPLICITIES...> impl_t;
          typedef eth::base::RefElList<MULTIPLICITIES::refElType...> refElList_t;


          /// This function should in fact be called, but Intel doesn't support this always...
          static constexpr eth::base::unsigned_t multiplicity2(eth::base::RefElType type) {
            return getMult<MULTIPLICITIES...>(type);
          }

          /// Unfortunately Intel Compiler doesn't allows us to call multiplicity directly,
          /// so we have to go over this helper...
          template<eth::base::RefElType TYPE>
          struct multiplicity {
            static const eth::base::unsigned_t val = multiplicity2(TYPE);
          };

        private:
          template<class MULT, class MULT2, class... MULTIPLICITIES_>
          static constexpr eth::base::unsigned_t getMult(eth::base::RefElType type) {
            return MULT::refElType==type ? MULT::multiplicity : getMult<MULT2,MULTIPLICITIES_...>(type);
          }

          template<class MULT>
          static constexpr eth::base::unsigned_t getMult(eth::base::RefElType type) {
            return MULT::refElType==type ? MULT::multiplicity : throw std::out_of_range("The dataset does not attach data to this refElType.");
          }

        }; // struct GridDataSetTraits
      } // namespace detail
    } // namespace utils
  } // namespace grids
} // namespace eth



#include "i_grid_data_set.hpp"

namespace eth {
  namespace grids {
    namespace utils {
   
      /**  
       *  \brief An implementation of IGridDataSet which updates its values when  
       *         the grid changes.
       *  \tparam DATA_TYPE         The data type of the data entries attached to entities of the grid view
       *  \tparam GRID_VIEW_FACTORY A class which derives from `IGridDataSet` interface...
       *  \tparam MULTIPLICITIES... A sequence of MultiplicityPair structs defining the reference elements on
       *                            which data entries shall be attached, and the number of such data entries
       *                            for each reference elements.
       *                            
       * \note There exist many convenient typedefs for the GridDataSet which
       *       have a much more concise (but less general) syntax to create
       *       a gridDataSet. Usually it's easier to use them.
       *       
       * \note When the grid changes (e.g. through refinement) this GridDataSet
       *       class will automatically adjust its internal structure to fit to
       *       the new grid (this happens automatically, the user must \e not
       *       call any function). The implementation makes sure that all entities
       *       which have the same id (`eth::grid::IdSet::id()` ) before and after
       *       the change will have the same data attached to them after the change.
       *       If entities are removed, their data is lost and if entities
       *       are added they are either initialized to a default value or
       *       remain in an uninitialized state (depending on which constructor
       *       is used to create the dataset).
       *       In particular this means that the user must \e manually attach 
       *       custom data to newly created entities (new elements can easily 
       *       be identified through eth::grid::Entity< GRID_TRAITS, 0 >::isNew()
       *       and from this one can deduce whether its subentities are new).
       *       
       */
      template<class DATA_TYPE, class GRID_VIEW_FACTORY, class... MULTIPLICITIES>
      class GridDataSet : public IGridDataSet<detail::GridDataSetTraits<DATA_TYPE,GRID_VIEW_FACTORY,MULTIPLICITIES...>> {
        // forward declarations:
      private:
        template<eth::base::RefElType TYPE, eth::base::unsigned_t MULTIPLICITY>
        class GDSHelper;
      public:
        /// The grid data set traits class
        typedef detail::GridDataSetTraits<DATA_TYPE,GRID_VIEW_FACTORY,MULTIPLICITIES...> traits_t;
        /// The data type of the data entries attached to the entities of the grid view
        typedef DATA_TYPE dataType_t;
        /// A model of \ref GridTraitsDoc defining the type of the grid view to whose entities data shall be attached
        typedef typename traits_t::viewTraits_t viewTraits_t;
        /// A model of \ref ViewTraitsDoc corresponding to viewTraits_t
        typedef typename viewTraits_t::gridTraits_t gridTraits_t;
        /// The base class type (CRTP)
        typedef IGridDataSet<traits_t> base_t;
        /// The grid view factory type required by the constructor to be able to re-create grid views even after mesh updates
        typedef GRID_VIEW_FACTORY gridViewFactory_t;
        /// The type of index set returned by the grid view
        typedef eth::grid::IndexSet<gridTraits_t,typename viewTraits_t::indexSet_t> indexSet_t;
        /// The composite type in which the data are stored per reference element
        typedef std::tuple<GDSHelper<MULTIPLICITIES::refElType, MULTIPLICITIES::multiplicity>...> tuple_t;
        typedef typename gridTraits_t::size_type size_type;
        /// A list of all reference element types to which this grid data set attaches data.
        typedef typename eth::base::RefElList<MULTIPLICITIES::refElType...> refElList_t;

        /** \brief Constructor \e without default value
         *  \param gridViewFactory GridViewFactory object which provides a nullary function to obtain a specific grid view of a specified grid
         *                         When the grid changes, new entries are not initialized.
         *  \note Upon creation of the GridDataSet, the values are not initialized and if new entities are added, they
         *        are also not set to a default value.
         */
        GridDataSet( gridViewFactory_t gridViewFactory ) :
          base_t(),
          gridViewFactory_( std::move(gridViewFactory) ),
          gridView_( std::move(gridViewFactory_.getView()) ),
          indexSet_( &gridView_.indexSet() ),
          initializeOnResize_( false )
        {
          resizeGDSHelper<MULTIPLICITIES...>();
        }

        /** \brief Constructor with default value
         *  \param gridViewFactory  Object which provides a nullary function to obtain a specific grid view of a specified grid
         *  \param defaultValue     Default value to which data is initialized and to which new data entries are initialized if the grid changes.
         */
        GridDataSet( gridViewFactory_t gridViewFactory, 
                     const dataType_t& defaultValue ) :
          base_t(),
          gridViewFactory_( std::move(gridViewFactory) ),
          gridView_( std::move(gridViewFactory_.getView()) ),
          indexSet_( &gridView_.indexSet() ),
          initializeOnResize_( true ),
          defaultValue_( defaultValue )
        {
          resizeGDSHelper<MULTIPLICITIES...>(defaultValue_);
        }

        /// destructor
        ~GridDataSet() {}

        /// see IGridDataSet::data(const eth::grid::Entity<gridTraits_t,CODIM>& e, size_type i)
        template<int CODIM>
        dataType_t& data(const eth::grid::Entity<gridTraits_t,CODIM>& e, size_type i) {
          static_assert(doWeSupportDim<MULTIPLICITIES...>(gridTraits_t::dimMesh-CODIM),
            "GridDataSet does not attach data to entities of this codim.");
          ETH_ASSERT_MSG ( doWeSupportRefElType<MULTIPLICITIES...>( e.refElType() ), 
            (std::string("RefElType ") + eth::base::getRefElName( e.refElType() )
            + std::string(" is not supported.") ).c_str());
          ETH_ASSERT_MSG( i < traits_t::multiplicity2( e.refElType() ),
            ( std::string("RefElType ") + eth::base::getRefElName( e.refElType() ) +
            std::string(" has only ") + 
            boost::lexical_cast<std::string>(traits_t::multiplicity2( e.refElType() ) ) +
            std::string(" data entries attached, but you asked for entry ") +
            boost::lexical_cast<std::string>(i) ).c_str() );
          static constexpr int dim = gridTraits_t::dimMesh-CODIM;
          return eth::base::ReferenceElements::applyUniversal<DataEntryRwFunctor,typename refElList_t::template onlyDim_t<dim>,GridDataSet&,const eth::grid::Entity<gridTraits_t,CODIM>&,size_type>(e.refElType(),*this,e,i);
        }

        /// see documentation of IGridDataSet
        template<int CODIM>
        dataType_t& subData(const eth::grid::Entity<gridTraits_t,0>& e, size_type subIndex, size_type i) {
          static_assert(doWeSupportDim<MULTIPLICITIES...>(gridTraits_t::dimMesh-CODIM),
            "GridDataSet does not attach data to entities of this codim.");
          ETH_ASSERT_MSG ( doWeSupportRefElType<MULTIPLICITIES...>( eth::base::ReferenceElements::getSubEntityType( e.refElType(), CODIM, subIndex ) ), 
            (std::string("RefElType ") + eth::base::getRefElName( eth::base::ReferenceElements::getSubEntityType( e.refElType(), CODIM, subIndex ) )
            + std::string(" is not supported.") ).c_str());
          ETH_ASSERT_MSG( i < traits_t::multiplicity2( eth::base::ReferenceElements::getSubEntityType( e.refElType(), CODIM, subIndex ) ),
            ( std::string("RefElType ") + eth::base::getRefElName( eth::base::ReferenceElements::getSubEntityType( e.refElType(), CODIM, subIndex ) ) +
            std::string(" has only ") + 
            boost::lexical_cast<std::string>(traits_t::multiplicity2( eth::base::ReferenceElements::getSubEntityType( e.refElType(), CODIM, subIndex ) ) ) +
            std::string(" data entries attached, but you asked for entry ") +
            boost::lexical_cast<std::string>(i) ).c_str() );
          static constexpr int dim = gridTraits_t::dimMesh-CODIM;
          return eth::base::ReferenceElements::applyUniversal<SubDataEntryRwFunctor,typename refElList_t::template onlyDim_t<dim>,GridDataSet&,const eth::grid::Entity<gridTraits_t,0>&,size_type,size_type>(eth::base::ReferenceElements::getSubEntityType( e.refElType(), CODIM, subIndex ),*this,e,subIndex,i);
        }

        /// see IGridDataSet::data(const eth::grid::Entity<gridTraits_t,CODIM>& e, size_type i) const
        template<int CODIM>
        const dataType_t& data(const eth::grid::Entity<gridTraits_t,CODIM>& e, size_type i) const {
          static_assert(doWeSupportDim<MULTIPLICITIES...>(gridTraits_t::dimMesh-CODIM),
            "GridDataSet does not attach data to entities of this codim.");
          ETH_ASSERT_MSG ( doWeSupportRefElType<MULTIPLICITIES...>( e.refElType() ), 
            (std::string("RefElType ") + eth::base::getRefElName( e.refElType() )
            + std::string(" is not supported.") ).c_str());
          ETH_ASSERT_MSG( i < traits_t::multiplicity2( e.refElType() ),
            ( std::string("RefElType ") + eth::base::getRefElName( e.refElType() ) +
            std::string(" has only ") + 
            boost::lexical_cast<std::string>(traits_t::multiplicity2( e.refElType() ) ) +
            std::string(" data entries attached, but you asked for entry ") +
            boost::lexical_cast<std::string>(i) ).c_str() );
          static const int dim = gridTraits_t::dimMesh-CODIM;
          return eth::base::ReferenceElements::applyUniversal<DataEntryRoFunctor,typename refElList_t::template onlyDim_t<dim>,const GridDataSet&,const eth::grid::Entity<gridTraits_t,CODIM>&,size_type>(e.refElType(),*this,e,i);
        }

        /// see documentation of IGridDataSet
        template<int CODIM>
        const dataType_t& subData(const eth::grid::Entity<gridTraits_t,0>& e, size_type subIndex, size_type i) const {
          static_assert(doWeSupportDim<MULTIPLICITIES...>(gridTraits_t::dimMesh-CODIM),
            "GridDataSet does not attach data to entities of this codim.");
          ETH_ASSERT_MSG ( doWeSupportRefElType<MULTIPLICITIES...>( eth::base::ReferenceElements::getSubEntityType( e.refElType(), CODIM, subIndex ) ), 
            (std::string("RefElType ") + eth::base::getRefElName( eth::base::ReferenceElements::getSubEntityType( e.refElType(), CODIM, subIndex ) )
            + std::string(" is not supported.") ).c_str());
          ETH_ASSERT_MSG( i < traits_t::multiplicity2( eth::base::ReferenceElements::getSubEntityType( e.refElType(), CODIM, subIndex ) ),
            ( std::string("RefElType ") + eth::base::getRefElName( eth::base::ReferenceElements::getSubEntityType( e.refElType(), CODIM, subIndex ) ) +
            std::string(" has only ") + 
            boost::lexical_cast<std::string>(traits_t::multiplicity2( eth::base::ReferenceElements::getSubEntityType( e.refElType(), CODIM, subIndex ) ) ) +
            std::string(" data entries attached, but you asked for entry ") +
            boost::lexical_cast<std::string>(i) ).c_str() );
          static const int dim = gridTraits_t::dimMesh-CODIM;
          return eth::base::ReferenceElements::applyUniversal<SubDataEntryRoFunctor,typename refElList_t::template onlyDim_t<dim>,const GridDataSet&,const eth::grid::Entity<gridTraits_t,0>&,size_type,size_type>(eth::base::ReferenceElements::getSubEntityType( e.refElType(), CODIM, subIndex ),*this,e,subIndex,i);
        }

        /// see documentation in IGridDataSet
        template<eth::base::RefElType TYPE>
        std::array<dataType_t,traits_t::template multiplicity<TYPE>::val>& data
          (const eth::grid::Entity<gridTraits_t,gridTraits_t::dimMesh-eth::base::ReferenceElement<TYPE>::dimension>& e) {
            ETH_ASSERT_MSG ( doWeSupportRefElType<MULTIPLICITIES...>( e.refElType() ), 
            (std::string("RefElType ") + eth::base::getRefElName( e.refElType() )
            + std::string(" is not supported.") ).c_str());
            return static_cast<GDSHelper<TYPE,traits_t::template multiplicity<TYPE>::val>*>(getGDSHelper<TYPE,MULTIPLICITIES...>())->data(e);
        }

        /// see documentation in IGridDataSet
        template<eth::base::RefElType TYPE>
        std::array<dataType_t,traits_t::template multiplicity<TYPE>::val>& subData
          (const eth::grid::Entity<gridTraits_t,0>& e, size_type subIndex) {
            ETH_ASSERT_MSG ( doWeSupportRefElType<MULTIPLICITIES...>( eth::base::ReferenceElements::getSubEntityType( e.refElType(), gridTraits_t::dimMesh-eth::base::ReferenceElement<TYPE>::dimension, subIndex ) ), 
              (std::string("RefElType ") + eth::base::getRefElName( eth::base::ReferenceElements::getSubEntityType( e.refElType(), gridTraits_t::dimMesh-eth::base::ReferenceElement<TYPE>::dimension, subIndex ) )
              + std::string(" is not supported.") ).c_str());
            return static_cast<GDSHelper<TYPE,traits_t::template multiplicity<TYPE>::val>*>(getGDSHelper<TYPE,MULTIPLICITIES...>())->subData(e,subIndex);
        }

        /// see documentation in IGridDataSet
        template<eth::base::RefElType TYPE>
        const std::array<dataType_t,traits_t::template multiplicity<TYPE>::val>& data
          (const eth::grid::Entity<gridTraits_t,gridTraits_t::dimMesh-eth::base::ReferenceElement<TYPE>::dimension>& e) const {
          ETH_ASSERT_MSG ( doWeSupportRefElType<MULTIPLICITIES...>( e.refElType() ), 
            (std::string("RefElType ") + eth::base::getRefElName( e.refElType() )
            + std::string(" is not supported.") ).c_str());
            return static_cast<const GDSHelper<TYPE,traits_t::template multiplicity<TYPE>::val>*>(getGDSHelper<TYPE,MULTIPLICITIES...>())->data(e);
        }

        /// see documentation in IGridDataSet
        template<eth::base::RefElType TYPE>
        const std::array<dataType_t,traits_t::template multiplicity<TYPE>::val>& subData
          (const eth::grid::Entity<gridTraits_t,0>& e, size_type subIndex) const {
            ETH_ASSERT_MSG ( doWeSupportRefElType<MULTIPLICITIES...>( eth::base::ReferenceElements::getSubEntityType( e.refElType(), gridTraits_t::dimMesh-eth::base::ReferenceElement<TYPE>::dimension, subIndex ) ), 
              (std::string("RefElType ") + eth::base::getRefElName( eth::base::ReferenceElements::getSubEntityType( e.refElType(), gridTraits_t::dimMesh-eth::base::ReferenceElement<TYPE>::dimension, subIndex ) )
              + std::string(" is not supported.") ).c_str());
            return static_cast<const GDSHelper<TYPE,traits_t::template multiplicity<TYPE>::val>*>(getGDSHelper<TYPE,MULTIPLICITIES...>())->subData(e,subIndex);
        }


      private:

        /// Functor for read/write access to data entry
        template<eth::base::RefElType TYPE>
        struct DataEntryRwFunctor {
          typedef dataType_t& returnType_t;
          static returnType_t invoke(GridDataSet& ds,const eth::grid::Entity<gridTraits_t,gridTraits_t::dimMesh-eth::base::ReferenceElement<TYPE>::dimension>& e, size_type i) {
            return static_cast<GDSHelper<TYPE,traits_t::template multiplicity<TYPE>::val>*>(ds.getGDSHelper<TYPE,MULTIPLICITIES...>())->data(e,i);
          }
        };

        /// Functor for read/write access to subData entry
        template<eth::base::RefElType TYPE>
        struct SubDataEntryRwFunctor {
          typedef dataType_t& returnType_t;
          static returnType_t invoke(GridDataSet& ds,const eth::grid::Entity<gridTraits_t,0>& e, size_type subIndex, size_type i) {
            return static_cast<GDSHelper<TYPE,traits_t::template multiplicity<TYPE>::val>*>(ds.getGDSHelper<TYPE,MULTIPLICITIES...>())->subData(e,subIndex,i);
          }
        };

        /// Functor for read-only access to data entry
        template<eth::base::RefElType TYPE>
        struct DataEntryRoFunctor {
          typedef const dataType_t& returnType_t;
          static returnType_t invoke(const GridDataSet& ds,const eth::grid::Entity<gridTraits_t,gridTraits_t::dimMesh-eth::base::ReferenceElement<TYPE>::dimension>& e, size_type i) {
            return static_cast<const GDSHelper<TYPE,traits_t::template multiplicity<TYPE>::val>*>(ds.getGDSHelper<TYPE,MULTIPLICITIES...>())->data(e,i);
          }
        };

        /// Functor for read-only access to subData entry
        template<eth::base::RefElType TYPE>
        struct SubDataEntryRoFunctor {
          typedef const dataType_t& returnType_t;
          static returnType_t invoke(const GridDataSet& ds,const eth::grid::Entity<gridTraits_t,0>& e, size_type subIndex, size_type i) {
            return static_cast<const GDSHelper<TYPE,traits_t::template multiplicity<TYPE>::val>*>(ds.getGDSHelper<TYPE,MULTIPLICITIES...>())->subData(e,subIndex,i);
          }
        };

        /// non-template base class for template class GDSHelper to allow as return type if needed
        class IGDSHelper {};

        /** \brief class to store a MULTIPLICITY of data entries for each entity with reference element of type TYPE
         *  \tparam MULTIPLICITY the number of data entries per entity with reference element of type TYPE
         *  \tparam TYPE         the type of reference element of the entries for which data are to be stored
         */
        template<eth::base::RefElType TYPE, eth::base::unsigned_t MULTIPLICITY>
        class GDSHelper : public IGDSHelper {
        public:
          /// the type of reference element of the entries for which data are to be stored
          static const eth::base::RefElType refElType = TYPE;
          /// the codimension of the entities corresponding to the reference element of type TYPE
          static const int codim = gridTraits_t::dimMesh-eth::base::ReferenceElement<TYPE>::dimension;
          /// default constructor, required for use in std::tuple
          GDSHelper() :
            indexSet_( nullptr )
          {}

          /// resize data vector, e.g. after grid change
          void resize( const indexSet_t* indexSet ) {
            indexSet_ = indexSet;
            data_.resize( indexSet_->size(TYPE) );
          }

          /// resize data vector with initialization, e.g. after grid change
          void resize( const indexSet_t* indexSet,
                       const dataType_t& defaultValue ) {
            indexSet_ = indexSet;
            std::array<dataType_t,MULTIPLICITY> defaultValueArray;
            for(eth::base::unsigned_t i=0; i<MULTIPLICITY; ++i) {
              defaultValueArray[i] = defaultValue;
            }
            data_.resize( indexSet_->size(TYPE), defaultValueArray );
          }

          /// see IGridDataSet::data(const eth::grid::Entity<gridTraits_t,CODIM>& e, size_type i)
          dataType_t&
          data( const eth::grid::Entity<gridTraits_t,codim>& e,
                size_type i ) {
            return data_[indexSet_->index( e )][i];
          }

          /// see IGridDataSet::subData(const eth::grid::Entity<gridTraits_t,0>& e, size_type subIndex, size_type i)
          dataType_t&
          subData( const eth::grid::Entity<gridTraits_t,0>& e,
                   size_type subIndex,
                   size_type i ) {
            return data_[ indexSet_->template subIndex<codim>( e, subIndex ) ][i];
          }

          /// see IGridDataSet::data(const eth::grid::Entity<gridTraits_t,CODIM>& e, size_type i) const
          const dataType_t& 
          data( const eth::grid::Entity<gridTraits_t,codim>& e, 
                size_type i ) const {
            return data_[indexSet_->index( e )][i];
          }

          /// see IGridDataSet::subData(const eth::grid::Entity<gridTraits_t,0>& e, size_type subIndex, size_type i) const
          const dataType_t&
          subData( const eth::grid::Entity<gridTraits_t,0>& e,
                   size_type subIndex,
                   size_type i ) const {
            return data_[ indexSet_->template subIndex<codim>( e, subIndex ) ][i];
          }

          /// see documentation in IGridDataSet
          std::array<dataType_t,MULTIPLICITY>&
          data( const eth::grid::Entity<gridTraits_t,codim>& e ) {
            return data_[indexSet_->index( e )];
          }

          /// see documentation in IGridDataSet
          std::array<dataType_t,MULTIPLICITY>&
          subData( const eth::grid::Entity<gridTraits_t,0>& e,
                   size_type subIndex ) {
            return data_[indexSet_->template subIndex<codim>( e, subIndex )];
          }

          /// see documentation in IGridDataSet
          const std::array<dataType_t,MULTIPLICITY>&
          data( const eth::grid::Entity<gridTraits_t,codim>& e) const {
            return data_[indexSet_->index( e )];
          }

          /// see documentation in IGridDataSet
          const std::array<dataType_t,MULTIPLICITY>&
          subData( const eth::grid::Entity<gridTraits_t,0>& e,
                   size_type subIndex ) const {
            return data_[indexSet_->template subIndex<codim>( e, subIndex )];
          }

        private:
          const indexSet_t* indexSet_;
          std::vector<std::array<dataType_t,MULTIPLICITY>> data_;
        }; // class GridDataSet::GDSHelper

        template<class PAIR_1,class PAIR_2,class... PAIRS_>
        void resizeGDSHelper() {
          std::get<sizeof...(PAIRS_)+1>(gdsHelpers_).resize(indexSet_);
          resizeGDSHelper<PAIR_2,PAIRS_...>();
        }

        template<class PAIR>
        void resizeGDSHelper() {
          std::get<0>(gdsHelpers_).resize(indexSet_);
        }

        template<class PAIR_1,class PAIR_2,class... PAIRS_>
        void resizeGDSHelper( const dataType_t& defaultValue ) {
          std::get<sizeof...(PAIRS_)+1>(gdsHelpers_).resize(indexSet_, defaultValue);
          resizeGDSHelper<PAIR_2,PAIRS_...>(defaultValue);
        }

        template<class PAIR>
        void resizeGDSHelper( const dataType_t& defaultValue ) {
          std::get<0>(gdsHelpers_).resize(indexSet_, defaultValue);
        }

        template<class PAIR_1, class PAIR_2, class... PAIRS_>
        static constexpr bool doWeSupportRefElType( eth::base::RefElType type ) {
          return (type == PAIR_1::refElType) || doWeSupportRefElType<PAIR_2,PAIRS_...>( type );
        }

        template<class PAIR>
        static constexpr bool doWeSupportRefElType( eth::base::RefElType type ) {
          return type == PAIR::refElType;
        }

        /// Does the list of reference elements contain a reference element with the given dimension?
        template<class PAIR_1, class PAIR_2,class... PAIRS_>
        static constexpr bool doWeSupportDim(int dim) {
          return eth::base::ReferenceElement<PAIR_1::refElType>::dimension == dim ? 
            true : doWeSupportDim<PAIR_2,PAIRS_...>(dim);
        }
        template<class PAIR>
        static constexpr bool doWeSupportDim(int dim) {
          return eth::base::ReferenceElement<PAIR::refElType>::dimension == dim;
        }

        template<eth::base::RefElType TYPE, class PAIR_1, class PAIR_2, class...PAIRS_>
        IGDSHelper* getGDSHelper()  {
          return std::get<sizeof...(MULTIPLICITIES)-sizeof...(PAIRS_)-2>(gdsHelpers_).refElType==TYPE ?
            &std::get<sizeof...(MULTIPLICITIES)-sizeof...(PAIRS_)-2>(gdsHelpers_) :
          getGDSHelper<TYPE,PAIR_2, PAIRS_...>();
        }
        template<eth::base::RefElType TYPE, class PAIR>
        IGDSHelper* getGDSHelper() {
          return std::get<sizeof...(MULTIPLICITIES)-1>(gdsHelpers_).refElType==TYPE ? 
            &std::get<sizeof...(MULTIPLICITIES)-1>(gdsHelpers_) :
          throw std::logic_error("this should not be called.");
        }

        template<eth::base::RefElType TYPE, class PAIR_1, class PAIR_2, class...PAIRS_>
        const IGDSHelper* getGDSHelper() const  {
          return std::get<sizeof...(MULTIPLICITIES)-sizeof...(PAIRS_)-2>(gdsHelpers_).refElType==TYPE ?
            &std::get<sizeof...(MULTIPLICITIES)-sizeof...(PAIRS_)-2>(gdsHelpers_) :
          getGDSHelper<TYPE,PAIR_2, PAIRS_...>();
        }
        template<eth::base::RefElType TYPE, class PAIR>
        const IGDSHelper* getGDSHelper() const {
          return std::get<sizeof...(MULTIPLICITIES)-1>(gdsHelpers_).refElType==TYPE ? 
            &std::get<sizeof...(MULTIPLICITIES)-1>(gdsHelpers_) :
          throw std::logic_error("this should not be called.");
        }

        // DATA MEMBERS

        /// The data stored per reference element
        std::tuple<GDSHelper<MULTIPLICITIES::refElType, MULTIPLICITIES::multiplicity>...> gdsHelpers_;
        ///
        const gridViewFactory_t gridViewFactory_;
        eth::grid::GridView<viewTraits_t> gridView_;
        const indexSet_t* indexSet_;
        bool initializeOnResize_;
        dataType_t defaultValue_;

      }; // class GridDataSet

      namespace detail {
        /// THIS CODE DOES NOT COMPILE WITH INTEL, but it compiles with CLang and g++
        /// Todo: uncomment this code as soon as Intel supports it.
        /// Some helper structs which are needed for the typedefs:
        /*template<int CODIM,int MULTIPLICITY,class DATA_TYPE,class GRID_VIEW_FACTORY>
        struct singleCodimSingleMultiplicityGDS {
        private:
          template<eth::base::RefElType... TYPES_>
          struct helper {
            typedef GridDataSet<DATA_TYPE,GRID_VIEW_FACTORY,MultiplicityPair<TYPES_,MULTIPLICITY>...> gds_t;
          };
        public:
          typedef typename GRID_VIEW_FACTORY::gridTraits_t gridTraits_t;
          static const int dim = gridTraits_t::dimMesh-CODIM;
          typedef typename eth::base::AllRefElTypesList_t::onlyDim_t<dim>::template getParamPack<helper>::gds_t gds_t;
        };*/

        // INTEL WORKAROUND:
        template<int CODIM, int MULTIPLICITY, class DATA_TYPE, class GRID_VIEW_FACTORY>
        struct singleCodimSingleMultiplicityGDS {
        private:
          template<int DIM, class DUMMY> struct helper;
          template<class DUMMY> struct helper<0,DUMMY> {
            typedef GridDataSet<DATA_TYPE,GRID_VIEW_FACTORY,MultiplicityPair<eth::base::RefElType::POINT,MULTIPLICITY>> gds_t;
          };
          template<class DUMMY> struct helper<1,DUMMY> {
            typedef GridDataSet<DATA_TYPE,GRID_VIEW_FACTORY,MultiplicityPair<eth::base::RefElType::SEGMENT,MULTIPLICITY>> gds_t;
          };
          template<class DUMMY> struct helper<2,DUMMY> {
            typedef GridDataSet<DATA_TYPE,GRID_VIEW_FACTORY,MultiplicityPair<eth::base::RefElType::TRIA,MULTIPLICITY>,MultiplicityPair<eth::base::RefElType::QUAD,MULTIPLICITY>> gds_t;
          };
          template<class DUMMY> struct helper<3,DUMMY> {
            typedef GridDataSet<DATA_TYPE,GRID_VIEW_FACTORY,MultiplicityPair<eth::base::RefElType::TETRA,MULTIPLICITY>,MultiplicityPair<eth::base::RefElType::PYRAMID,MULTIPLICITY>,MultiplicityPair<eth::base::RefElType::PRISM,MULTIPLICITY>,MultiplicityPair<eth::base::RefElType::HEXA,MULTIPLICITY>> gds_t;
          };

        public:
          typedef typename helper<GRID_VIEW_FACTORY::gridTraits_t::dimMesh-CODIM,int>::gds_t gds_t;
        };
      }

      // Some convenient typedefs for GridDataSet:
      ////////////////////////////////////////////////////////////////////////// 
      
      /**
       * \brief A GridDataSet which attaches data to entities of a given codimension,
       *        Multiplicity =1
       */
      template<int CODIM, class DATA_TYPE, class GRID_VIEW_FACTORY >
      using SingleCodimGDS = typename detail::singleCodimSingleMultiplicityGDS<CODIM,1,DATA_TYPE,GRID_VIEW_FACTORY>::gds_t;

      /// The same as SingleCodimGDS, but here in terms of the interface class IGridDataSet,
      // useful to instantiate a shared_ptr with make_shared_crtp...
      template<int CODIM, class DATA_TYPE, class GRID_VIEW_FACTORY >
      using ISingleCodimGDS = typename detail::singleCodimSingleMultiplicityGDS<CODIM,1,DATA_TYPE,GRID_VIEW_FACTORY>::gds_t::base_t;

      

    } // namespace utils
  } // namespace grids
} // namespace eth

#endif // ETH_GRID_DATA_SET_HPP
