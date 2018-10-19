

#ifndef ETH_BOUNDARY_DATA_SET_HPP
#define ETH_BOUNDARY_DATA_SET_HPP

// own includes ----------------------------------------------------------------
#include <eth_base/ETH_ASSERT.hpp>

namespace eth {
  namespace grids {
    namespace utils {

      // Forward declarations:
      template<class DATA_TYPE, class GRID_TRAITS>
      class BoundaryDataSet;

      namespace detail {
        /// Traits class for the BoundaryDataSet class.
        template<class DATA_TYPE, class GRID_TRAITS>
        struct BoundaryDataSetTraits {
          typedef DATA_TYPE dataType_t;
          typedef GRID_TRAITS gridTraits_t;
          typedef BoundaryDataSet<dataType_t,gridTraits_t> impl_t;
        };
      } // namespace detail
    } // namespace utils
  } // namespace grids
} //namespace eth

#include "i_boundary_data_set.hpp"
#include <vector>

namespace eth {
  namespace grids {
    namespace utils {
      /**
       * \brief A (default) implementation of IBoundaryDataSet static interface
       *        which updates its attached data when the grid changes.
       * \tparam DATA_TYPE The type of data which is attached to boundary intersections of the grid.
       * \tparam GRID_TRAITS A type which models \ref GridTraitsDoc and specifies the underlying grid data structure.
       * 
       * \note When the grid changes (e.g. through refinement) this BoundaryDataSet
       *       class will automatically adjust its internal structure to fit to
       *       the new grid (this happens automatically, the user must \e not
       *       call any function). The implementation makes sure that all boundary intersections
       *       which have the same id (`eth::grid::IdSet::id()` ) before and after
       *       the change will have the same data attached to them after the change.
       *       If boundary intersections are removed, their data is lost and if intersections
       *       are added they are either initalized to a default value or
       *       remain in an uninitialized state (depending on which constructor
       *       is used to create the boundary dataset).
       *       In particular this means that the user must \e manually attach 
       *       custom data to newly created boundary intersections.
       */
      template<class DATA_TYPE, class GRID_TRAITS>
      class BoundaryDataSet : public IBoundaryDataSet<detail::BoundaryDataSetTraits<DATA_TYPE,GRID_TRAITS>> {
      public:
        typedef DATA_TYPE dataType_t;
        typedef GRID_TRAITS gridTraits_t;
        typedef detail::BoundaryDataSetTraits<dataType_t,gridTraits_t> traits_t;
        

        /**
         * \brief Default constructor without default value.
         *
         * \param [in,out]  grid  A shared pointer to the underlying grid
         *                  
         * When you call this constructor, the values attached to the boundary
         * intersections will remain in an uninitialized state and when new
         * boundary intersections are added to the grid they will also remain
         * in an invalid state.
         */
        BoundaryDataSet(std::shared_ptr<eth::grid::Grid<gridTraits_t>> grid) : 
          grid_(grid), 
          initializeOnResize_(false), 
          boundaryIndexSet_(grid->boundaryIndexSet()),
          data_(boundaryIndexSet_.size()){

        }

        /**
         * \brief Constructor which provides a default value with which
         *        new intersections are initialized.
         *
         * \param [in,out]  grid  A shared pointer to the underlying grid
         *                  
         * When you call this constructor, the values attached to the boundary
         * intersections will be initialized with the given value and if new
         * intersections are added later on they will also be initialized
         * with this default value.
         */
        BoundaryDataSet(std::shared_ptr<eth::grid::Grid<gridTraits_t>> grid, const dataType_t& defaultValue) :
          grid_(grid),
          initializeOnResize_(true),
          defaultValue_( defaultValue ),
          boundaryIndexSet_(grid->boundaryIndexSet()),
          data_(boundaryIndexSet_.size(),defaultValue_) {
        }

        /// destructor
        ~BoundaryDataSet() {}

        dataType_t& data(const eth::grid::Intersection<gridTraits_t>& i) {
          ETH_ASSERT_MSG(i.boundary(), "Data can only be stored with boundary intersections!!!");
          return data_[boundaryIndexSet_.index(i)];
        }

        const dataType_t& data(const eth::grid::Intersection<gridTraits_t>& i) const {
          ETH_ASSERT_MSG(i.boundary(), "Data can only be stored with boundary intersections!!!");
          return data_[boundaryIndexSet_.index(i)];
        }

      protected:
        /// shared pointer of the underlying grid...
        std::shared_ptr<eth::grid::Grid<gridTraits_t>> grid_;
        /// should we initialize data for newly added intersections?
        bool initializeOnResize_;
        /// if initializeOnResize_== true, we store here the actual value with
        // which new intersections are initialized.
        dataType_t defaultValue_;

        const eth::grid::BoundaryIndexSet<gridTraits_t>& boundaryIndexSet_;

        /// store data associated with every boundary intersection (through boundyIndexSet)
        std::vector<dataType_t> data_;
      };

    } // namespace utils
  } // namespace grids
} // namespace eth






#endif // ETH_BOUNDARY_DATA_SET_HPP
