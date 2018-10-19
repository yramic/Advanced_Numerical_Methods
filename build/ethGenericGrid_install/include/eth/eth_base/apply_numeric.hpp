/********************************************************************
    created:    2013/10/10
    created:    10:10:2013   17:20
    filename:   apply_numeric.hpp
    author:     Raffael Casagrande

    purpose:    A tool to convert a run-time int value into a compile
                time int value.
*********************************************************************/
//          Copyright Raffael Casagrande 2013 - 2015.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)


#ifndef ETH_BASE_APPLY_NUMERIC_HPP
#define ETH_BASE_APPLY_NUMERIC_HPP

#include <boost/lexical_cast.hpp>
#include "ETH_ASSERT.hpp"

namespace eth {
  namespace base {

    namespace detail {
      template<template<int> class FUNCTOR, int CURRENT, int NUMERIC_MIN, int NUMERIC_MAX, class... ARGS>
      auto applyNumericHelper(int a, ARGS... args) -> typename FUNCTOR<NUMERIC_MIN>::returnType_t {
        if(a == CURRENT)
          return FUNCTOR<CURRENT>::invoke(args...);
        else if(a<= NUMERIC_MAX)
          return applyNumericHelper<FUNCTOR,CURRENT+1,NUMERIC_MIN,NUMERIC_MAX,ARGS...>(a,args...);
        else
          ETH_VERIFY_MSG(false,(std::string("cannot call functor because a = ") +
          boost::lexical_cast<std::string>(a) + std::string(" is not in the range ") + 
          boost::lexical_cast<std::string>(NUMERIC_MIN) + std::string(" - ") +
          boost::lexical_cast<std::string>(NUMERIC_MAX)).c_str());
      }
    } // end namespace anonymous


    /**
     * \brief Helper template function which allows a runtime integer parameter
     *        to be converted into a compile time parameter through a "switch statement".
     * \tparam  FUNCTOR A functor which should be called (and to which the compile-time
     *                  argument will be provided).
     * \tparam NUMERIC_MIN  Lower bound in which the run-time value a must lie.
     * \tparam NUMERIC_MAX  Upper bound in which the run-time value a must lie.
     * \tparam ARGS         The types of the additional arguments which are passed
     *                      to the functor.
     * \param a     The runtime int parameter which should be converted.
     * \param args  additional arguments which are passed to the functor.
     * \return  The value returned by FUNCTOR<a>::invoke(args...)
     * 
     * \note The FUNCTOR must provide the public type `returnType_t`
     *       and a \e static function
     *       `returnType_t invoke(ARGS...)`
     *       
     * Conceptually the call to this function is the same as:
     * \code
     *      switch(a) {
     *        case NUMERIC_MIN:
     *          return FUNCTOR<NUMERIC_MIN>::invoke(args...);
     *        case NUMERIC_MIN+1:
     *          return FUNCTOR<NUMERIC_MIN+1>::invoke(args...);
     *        .
     *        .
     *        .
     *        case NUMERIC_MAX:
     *          return FUNCTOR<NUMERIC_MAX>::invoke(args...);
     *      }
     * \endcode
     */
    template<template<int> class FUNCTOR, int NUMERIC_MIN, int NUMERIC_MAX, class... ARGS>
    auto applyNumeric(int a, ARGS... args) -> typename FUNCTOR<NUMERIC_MIN>::returnType_t {
      return detail::applyNumericHelper<FUNCTOR,NUMERIC_MIN,NUMERIC_MIN,NUMERIC_MAX, ARGS...>(a,args...);
    }

    
  } // end namespace base
} // end namespace eth

#endif // ETH_BASE_APPLY_NUMERIC_HPP
