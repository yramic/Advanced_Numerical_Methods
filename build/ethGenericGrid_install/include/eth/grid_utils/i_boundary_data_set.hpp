#ifndef ETH_I_BOUNDARY_DATA_SET_HPP
#define ETH_I_BOUNDARY_DATA_SET_HPP

#include "interface.hpp"

namespace eth {
  namespace grids {
    namespace utils {
      
      /**
       * \brief A CRTP interface class which defines how a user can attach custom
       *        data to boundary intersections.
       * \tparam BDS_TRAITS A class which defines implementation specific types
       *                    (see below for more details.)
       * 
       * A IBoundaryDataSet::data allows the user to attach any kind of data
       * (type: dataType_t) to boundary intersections 
       * (i.e. `eth::grid::Intersection::boundary()==true`)
       * 
       * \note A IBoundaryDataSet is neither copyable nor moveable. If you want
       *       to create a deep copy of a BoundaryDataSet you must cast it into
       *       the actual implementation (CRTP pattern) and call then  it's copy
       *       constructor.
       * 
       *
       * <H3> Note for Implementors </H3>
       * The BDS_TRAITS class must provide the following public types:
       * type                      | meaning
       * --------------------------|----------------------------------------------------
       * impl_t                    | The actual class which implements the IBoundaryDataSet interface (CRTP pattern)
       * dataType_t                | The type of data which is attached to the entities.
       * gridTraits_t              | A class which models \ref GridTraitsDoc and specifies the underlying grid data structure.
       */
      template<class BDS_TRAITS>
      class IBoundaryDataSet {
      public:
        typedef BDS_TRAITS traits_t;
        /// underlying CRTP implementation
        typedef typename traits_t::impl_t impl_t;
        /// The type of data which is attached to the entities.
        typedef typename traits_t::dataType_t dataType_t;
        /// The gridTraits class to which this BoundaryDataSet refers.
        typedef typename traits_t::gridTraits_t gridTraits_t;
        typedef eth::base::unsigned_t size_type;

        typedef typename gridTraits_t::intersection_t intersection_t;

        /**
         * \brief Retrieve (and possibly modify) data attached to \e boundary intersection i.
         * \param i Intersection object to which data is attached.
         * \warning if `i.boundary()==false` an assert statement will break the
         *          execution of the code.
         */
        dataType_t& data(const eth::grid::Intersection<gridTraits_t>& i) {
          CHECK_CRTP_RETURN(asImp().data(i));
        }

        /// const version of data(const eth::grid::Intersection<gridTraits_t>& i)
        const dataType_t& data(const eth::grid::Intersection<gridTraits_t>& i) const {
          CHECK_CRTP_RETURN(asImp().data(i));
        }

      protected:
        /// hide default Constructor (CRTP)
        IBoundaryDataSet() {}
        /// hide copy constructor (CRTP)
        IBoundaryDataSet(const IBoundaryDataSet&) {}
        // hide destructor (CRTP)
        ~IBoundaryDataSet() {}

      private:
        /// CRTP trick
        impl_t& asImp() { return static_cast<impl_t&> (*this);}
        /// CRTP trick
        const impl_t& asImp()const {return static_cast<const impl_t&>(*this);}
      };
    }
  }
}


#endif // ETH_I_BOUNDARY_DATA_SET_HPP
