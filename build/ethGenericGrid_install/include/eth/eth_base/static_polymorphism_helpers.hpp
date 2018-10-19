/********************************************************************
    created:    2013/10/21
    created:    21:10:2013   12:09
    filename:   static_polymorphism_helpers.hpp
    author:     Raffael Casagrande, Christoph Winkelmann

    purpose:    Helpers for static polymorphism
*********************************************************************/
//          Copyright Raffael Casagrande 2013 - 2015.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#ifndef _HPP_772FB358_5A4E_47DB_8504_2BD552F9C015
#define _HPP_772FB358_5A4E_47DB_8504_2BD552F9C015

namespace eth {
  namespace base {

    /**
     * \brief Test if a given class implements a given templated interface
     * \tparam IMPL      The interface class to test implementation.
     * \tparam INTERFACE The interface templatized with a traits class.
     */
    template<class IMPL,template<class TRAITS> class INTERFACE>
    class InterfaceTester {
      template<class TRAITS>
      static TRAITS  getTraits(INTERFACE<TRAITS>* /* a */) {
        return TRAITS();
      }

      static bool getTraits(void* /* a */) {
        return false;
      }

    public:
      /// true if IMPL is / implements INTERFACE
      static const bool isInterface = !std::is_same<decltype(getTraits((IMPL*)0)),bool>::value;

      /// the traits class with which IMPL implements INTERFACE
      typedef decltype(getTraits((IMPL*)0)) traits_t;
    };
    

  } // end namespace base
} // end namespace eth

#endif // _HPP_772FB358_5A4E_47DB_8504_2BD552F9C015
