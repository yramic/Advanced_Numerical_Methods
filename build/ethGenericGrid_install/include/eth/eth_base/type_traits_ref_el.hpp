/********************************************************************
  created:	2013/11/15
  created:	23:7:2013   15:04
  filename: type_traits_ref_el.hpp
  author:		Lars Kielhorn
  
  purpose:	Helper traits to define sub-entities on reference elements
*********************************************************************/

#ifndef ETH_GRID_TYPE_TRAITS_REF_EL_HPP
#define ETH_GRID_TYPE_TRAITS_REF_EL_HPP

// own include -----------------------------------------------------------------
#include "numeric_types.hpp"
#include "integer_list.hpp"
#include "ETH_ASSERT.hpp"
#include "enum_ref_el_types.hpp"

// boost includes --------------------------------------------------------------
#include <boost/lexical_cast.hpp>


namespace eth {
  namespace base {
    namespace detail {
      // Helper traits class to determine dimension of Reference Elements...
      // intel compiler support explicit specialization but g++, clang 
      // deal only with partial specialization of member classes. Hence, we add
      // the artificial template parameter DUMMY. Not nice! :-(
      template<RefElType TYPE2,int DUMMY=0>
      struct typeTraits_;


      //------------------------------------------------------------------------
      // POINT
      template<int DUMMY>
      struct typeTraits_<RefElType::POINT,DUMMY>
      {
        typedef integer_list<1,0,0,0,0,0,0,0> entity_list;
        static constexpr int dim = 0;
        template< int CODIM >
        static constexpr unsigned_t numSubEntities( ) {
          return entity_list::get<dim-CODIM>();
        }
        static unsigned_t numSubEntities(int codim) {
          return ( codim == 0 ? 
                   numSubEntities<0>() : 
                   static_cast< unsigned_t >( -1 ) );
        }
      };

      //------------------------------------------------------------------------
      // SEGMENT
      template<int DUMMY>
      struct typeTraits_<RefElType::SEGMENT,DUMMY>
      {
        typedef integer_list<2,1,0,0,0,0,0,0> entity_list;
        static constexpr int dim = 1;
        template< int CODIM >
        static constexpr unsigned_t numSubEntities( ) {
          return entity_list::get<dim-CODIM >( );
        }
        static unsigned_t numSubEntities(int codim) {
          switch( codim ) {
          case 0: return numSubEntities<0>( );
          case 1: return numSubEntities<1>( );
          default: return static_cast< unsigned_t >(-1);
          }
        }
      };


      //------------------------------------------------------------------------
      // TRIA
      template<int DUMMY>
      struct typeTraits_<RefElType::TRIA,DUMMY>
      { 
        typedef integer_list<3,3,1,0,0,0,0,0> entity_list;
        static constexpr int dim = 2;
        template< int CODIM >
        static constexpr unsigned_t numSubEntities( ) {
          return entity_list::get<dim-CODIM>( );
        }
        static unsigned_t numSubEntities(int codim) {
          switch( codim ) {
          case 0: return numSubEntities<0>( );
          case 1: return numSubEntities<1>( ); 
          case 2: return numSubEntities<2>( );
          default: return static_cast< unsigned_t >( -1 );
          }
        }
      };
    

      //------------------------------------------------------------------------
      // QUAD
      template<int DUMMY>
      struct typeTraits_<RefElType::QUAD,DUMMY>
      { 
        typedef integer_list<4,4,0,1,0,0,0,0> entity_list;
        static constexpr int dim = 2;
        template< int CODIM >
        static constexpr unsigned_t numSubEntities( ) {
          return entity_list::get<dim-CODIM2Index_<CODIM>::mapped_codim>();
        }
        static unsigned_t numSubEntities(int codim) {
          switch( codim ) {
          case 0: return numSubEntities<0>( );
          case 1: return numSubEntities<1>( ); 
          case 2: return numSubEntities<2>( );
          default: return static_cast< unsigned_t >( -1 );
          }
        }

      private:
        template< int CODIM, int DUMMY2 = 0 >
        struct CODIM2Index_ { static const int mapped_codim = CODIM; };

        template< int DUMMY2 >
        struct CODIM2Index_< 0, DUMMY2 > { static const int mapped_codim = -1; };
      };



      //------------------------------------------------------------------------
      // TETRA
      template<int DUMMY>
      struct typeTraits_<RefElType::TETRA,DUMMY>
      { 
        typedef integer_list<4,6,4,0,1,0,0,0> entity_list;
        static constexpr int dim = 3;
        template< int CODIM >
        static constexpr unsigned_t numSubEntities( ) {
          return entity_list::get<dim-CODIM2Index_<CODIM>::mapped_codim>();
        }
        static unsigned_t numSubEntities(int codim) {
          switch( codim ) {
          case 0: return numSubEntities<0>( ); break;
          case 1: return numSubEntities<1>( ); break; 
          case 2: return numSubEntities<2>( ); break;
          case 3: return numSubEntities<3>( ); break;
          default: return static_cast< unsigned_t >( -1 );
          }
        }

      private:
        template< int CODIM, int DUMMY2 = 0 >
        struct CODIM2Index_ { static const int mapped_codim = CODIM; };

        template< int DUMMY2 >
        struct CODIM2Index_< 0, DUMMY2 > { static const int mapped_codim = -1; };
      };



      //------------------------------------------------------------------------
      // HEXA
      template<int DUMMY>
      struct typeTraits_<RefElType::HEXA,DUMMY>
      { 
        typedef integer_list<8,12,0,6,0,1,0,0> entity_list;
        static constexpr int dim = 3;

        template< int CODIM >
        static constexpr unsigned_t numSubEntities( ) {
          return entity_list::get<dim-CODIM2Index_<CODIM>::mapped_codim>();
        }
        static unsigned_t numSubEntities(int codim) {
          switch( codim ) {
          case 0: return numSubEntities<0>( ); break;
          case 1: return numSubEntities<1>( ); break;
          case 2: return numSubEntities<2>( ); break;
          case 3: return numSubEntities<3>( ); break;
          default: return static_cast< unsigned_t >( -1 );
          }
        }

      private:
        template< int CODIM, int DUMMY2 = 0 >
        struct CODIM2Index_ { static const int mapped_codim = CODIM; };

        template< int DUMMY2 >
        struct CODIM2Index_< 0, DUMMY2 > { static const int mapped_codim = -2; };

        template< int DUMMY2 >
        struct CODIM2Index_< 1, DUMMY2 > { static const int mapped_codim =  0; };
      };



      //------------------------------------------------------------------------
      // PYRAMID
      template<int DUMMY>
      struct typeTraits_<RefElType::PYRAMID,DUMMY>
      { 
        typedef integer_list<5,8,4,1,0,0,0,1> entity_list;
        static constexpr int dim = 3;
        template< int CODIM >
        static constexpr unsigned_t numSubEntities( ) {
          return CODIM2Entity_< CODIM >::value;
        }
        static unsigned_t numSubEntities(int codim) {
          switch( codim ) {
          case 0: return numSubEntities<0>( ); break;
          case 1: return numSubEntities<1>( ); break; 
          case 2: return numSubEntities<2>( ); break;
          case 3: return numSubEntities<3>( ); break;
          default: return static_cast< unsigned_t >( -1 );
          }
        }

      private:
        template< int CODIM, int DUMMY2 = 0 >
        struct CODIM2Entity_ { static const int value = 0; };

        template< int DUMMY2 >
        struct CODIM2Entity_< 0, DUMMY2 > { 
          static const int value = entity_list::get< 7 >( );
        };

        template< int DUMMY2 >
        struct CODIM2Entity_< 1, DUMMY2 > { 
          static const int value = entity_list::get<2>() + entity_list::get<3>();
        };

        template< int DUMMY2 >
        struct CODIM2Entity_< 2, DUMMY2 > { 
          static const int value = entity_list::get<1>();
        };

        template< int DUMMY2 >
        struct CODIM2Entity_< 3, DUMMY2 > { 
          static const int value = entity_list::get<0>();
        };
      };



      //------------------------------------------------------------------------
      // PRISM
      template<int DUMMY>
      struct typeTraits_<RefElType::PRISM,DUMMY>
      { 
        typedef integer_list<6,9,2,3,0,0,1,0> entity_list;
        static constexpr int dim = 3;
        template< int CODIM >
        static constexpr unsigned_t numSubEntities( ) {
          return CODIM2Entity_< CODIM >::value;
        }
        static unsigned_t numSubEntities(int codim) {
          switch( codim ) {
          case 0: return numSubEntities<0>( ); break;
          case 1: return numSubEntities<1>( ); break; 
          case 2: return numSubEntities<2>( ); break;
          case 3: return numSubEntities<3>( ); break;
          default: return static_cast< unsigned_t >( -1 );
          }
        }

      private:
        template< int CODIM, int DUMMY2 = 0 >
        struct CODIM2Entity_ { static const int value = 0; };

        template< int DUMMY2 >
        struct CODIM2Entity_< 0, DUMMY2 > { 
          static const int value = entity_list::get< 6 >( );
        };

        template< int DUMMY2 >
        struct CODIM2Entity_< 1, DUMMY2 > { 
          static const int value = entity_list::get<2>() + entity_list::get<3>();
        };

        template< int DUMMY2 >
        struct CODIM2Entity_< 2, DUMMY2 > { 
          static const int value = entity_list::get<1>();
        };

        template< int DUMMY2 >
        struct CODIM2Entity_< 3, DUMMY2 > { 
          static const int value = entity_list::get<0>();
        };
      };


    } // end namespace detail
  } // end namespace base
} // end namespace eth

#endif // ETH_GRID_TYPE_TRAITS_REF_EL_HPP
