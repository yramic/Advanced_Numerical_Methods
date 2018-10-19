#ifndef ETH_GRID_VIEW_FACTORY_HPP
#define ETH_GRID_VIEW_FACTORY_HPP

#include "interface.hpp"
#include "eth_base.hpp"
#include "eth_base/static_polymorphism_helpers.hpp"

namespace eth {
  namespace grids {
    namespace utils {
      
      /** \brief Provide nullary function to obtain a specific grid view of a specified grid
       *  \tparam GRID         A type which derives from grid::Grid (e.g. grids::hybrid::HybridGrid)
       *  \tparam TYPE         the type of the GridView which is returned
       *  This Factory allows the user to specify upon construction the grid, the type of the grid view and 
       *  the level (if applicable). The resulting grid view factory can then be passed onto
       *  any instance that repeatedly needs such grid views, also after changes of the grid.
       *  \note This is the general template which is not working as such. One of the two specializations 
       *        GridViewFactory<GRID_TRAITS, eth::grid::GridViewTypes::LeafView> or 
       *        GridViewFactory<GRID_TRAITS, eth::grid::GridViewTypes::LevelView> needs to be used.
       */
      template<class GRID, eth::grid::GridViewTypes TYPE>
      class GridViewFactory;
      
      /// Implementation for a leafView...
      template<class GRID>
      class GridViewFactory<GRID, eth::grid::GridViewTypes::LeafView> {
      public:

        static_assert(base::InterfaceTester<GRID,eth::grid::Grid>::isInterface,
          "The template parameter GRID must be a class deriving from grid::Grid<...>");

        static_assert(std::is_same<GRID,typename GRID::impl_t>::value,
          "The template parameter GRID must contain the implementation class, not the interface class.");

        /// the grid traits
        typedef typename GRID::gridTraits_t gridTraits_t;
        /// the grid view traits
        typedef typename gridTraits_t::template viewTraits_t<eth::grid::GridViewTypes::LeafView> viewTraits_t;
        /// the grid implementation
        using gridImpl_t = typename gridTraits_t::gridImpl_t;
        
        /** \brief Constructor
         * \param grid The grid to obtain leaf views from
         */
        GridViewFactory( std::shared_ptr<eth::grid::Grid<gridTraits_t>> grid )
          : grid_(grid) {}
        /// Obtain a leaf view of the grid which was specified upon construction
        eth::grid::GridView<viewTraits_t> getView() const {
          return grid_->leafView();
        }
        /// return the grid implementation
        const gridImpl_t& impl() const
        { return static_cast<const gridImpl_t&>( *grid_ ); }
        /// return the grid implementation (non-const version)
        gridImpl_t& impl() { return static_cast<gridImpl_t&>( *grid_ ); }

      private:
        std::shared_ptr<eth::grid::Grid<gridTraits_t>> grid_;
      };

      /** \brief Provide nullary function to obtain a level view of a specified grid
       *  \tparam GRID_TRAITS  type which models \ref GridTraitsDoc
       */
      template<class GRID>
      class GridViewFactory<GRID, eth::grid::GridViewTypes::LevelView> {
      public:

        static_assert(base::InterfaceTester<GRID,eth::grid::Grid>::isInterface,
          "The template parameter GRID must be a class deriving from grid::Grid<...>");

        static_assert(std::is_same<GRID,typename GRID::impl_t>::value,
          "The template parameter GRID must contain the implementation class, not the interface class.");


        /// the grid traits
        typedef typename GRID::gridTraits_t gridTraits_t;
        /// the grid view traits
        typedef typename gridTraits_t::template viewTraits_t<eth::grid::GridViewTypes::LevelView> viewTraits_t;
        /// the grid implementation
        using gridImpl_t = typename gridTraits_t::gridImpl_t;

        typedef typename gridTraits_t::size_type size_type;
        /** \brief constructor
         * \param grid The grid to obtain level views from
         * \param level The level of the view which is to be returned
         */
        GridViewFactory( std::shared_ptr<eth::grid::Grid<gridTraits_t>> grid, size_type level )
          : grid_(grid), level_(level) {}
        /// Obtain a level view of the grid which was specified upon construction
        eth::grid::GridView<viewTraits_t> getView() const {
          return grid_->levelView( level_ );
        }
        /// return the grid implementation
        const gridImpl_t& impl() const
        { return static_cast<const gridImpl_t&>( *grid_ ); }
        /// return the grid implementation (non-const version)
        gridImpl_t& impl() { return static_cast<gridImpl_t&>( *grid_ ); }
      private:
        std::shared_ptr<eth::grid::Grid<gridTraits_t>> grid_;
        size_type level_;
      };

    }
  }
}

#endif // ETH_GRID_VIEW_FACTORY_HPP
