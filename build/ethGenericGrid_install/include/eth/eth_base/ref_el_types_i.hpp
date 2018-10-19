#ifndef ETH_BASE_REF_EL_TYPES_I_HPP
#define ETH_BASE_REF_EL_TYPES_I_HPP

// own includes ---------------------------------------------------------------
#include "ref_el_types.hpp"
#include "numeric_types.hpp"

//-----------------------------------------------------------------------------
namespace eth {
  namespace base {

    namespace /* anonymous */ {
      // General Typedefs
      //////////////////////////////////////////////////////////////////////////
      /// Type used to count nodes and other things
      typedef eth::base::unsigned_t					size_type;
      /// Type used to group node/edge/face numbers together
      typedef std::vector<size_type>				indexVec_t;
      /// Type to group RefElTypes together...
      typedef std::vector<eth::base::RefElType>	typeVec_t;

      typedef std::vector<indexVec_t>				iiVec_t;
      typedef std::vector<typeVec_t>				itVec_t;
    } // end namespace anonymous

    //--------------------------------------------------------------------------
    // POINT
    //--------------------------------------------------------------------------
#if !defined (__INTEL_COMPILER)
    template<>
#endif
    const itVec_t ReferenceElement<RefElType::POINT>::
    subEntityTypes_{ {RefElType::POINT} };

#if !defined (__INTEL_COMPILER)
    template<>
#endif
    const std::vector<iiVec_t> ReferenceElement<RefElType::POINT>::
    subEntityCorners_{ { { 0 } } };
    //#endif // NO_INITIALIZER_LIST_SUPPORT

#if !defined (__INTEL_COMPILER)
    template<>
#endif
    const fixedMatrix_t<0,1> 
    ReferenceElement<RefElType::POINT>::cornerCoord_((fixedMatrix_t<0,1>()));

    //--------------------------------------------------------------------------
    // SEGMENT
    //--------------------------------------------------------------------------
#if !defined (__INTEL_COMPILER)
    template<>
#endif
    const itVec_t ReferenceElement<RefElType::SEGMENT>::subEntityTypes_{
      {RefElType::SEGMENT},{RefElType::POINT,RefElType::POINT}};

#if !defined (__INTEL_COMPILER)
    template<>
#endif
    const std::vector<iiVec_t> ReferenceElement<RefElType::SEGMENT>::subEntityCorners_{
      { {0,1} }, { {0},{1} } };
    //#endif // NO_INITIALIZER_LIST_SUPPORT

#if !defined (__INTEL_COMPILER)
    template<>
#endif
    const fixedMatrix_t<1,2>
    ReferenceElement<RefElType::SEGMENT>::cornerCoord_((fixedMatrix_t<1,2>() << -1,1).finished());


    //--------------------------------------------------------------------------
    // TRIA
    //--------------------------------------------------------------------------
#if !defined (__INTEL_COMPILER)
    template<>
#endif
    const itVec_t ReferenceElement<RefElType::TRIA>::subEntityTypes_{
      {RefElType::TRIA},{RefElType::SEGMENT,RefElType::SEGMENT,RefElType::SEGMENT},{RefElType::POINT,RefElType::POINT,RefElType::POINT}};

#if !defined (__INTEL_COMPILER)
    template<>
#endif
    const std::vector<iiVec_t> ReferenceElement<RefElType::TRIA>::subEntityCorners_{
      {{0,1,2}}, { {0,1},{1,2},{2,0} },{ {0},{1},{2} } };

#if !defined (__INTEL_COMPILER)
    template<>
#endif
    const fixedMatrix_t<2,3> ReferenceElement<RefElType::TRIA>::cornerCoord_(
                                                                             (fixedMatrix_t<2,3>() << 
                                                                              0,1,1,
                                                                              0,0,1).finished());



    //--------------------------------------------------------------------------
    // QUAD
    //--------------------------------------------------------------------------
#if !defined (__INTEL_COMPILER)
    template<>
#endif
    const itVec_t ReferenceElement<RefElType::QUAD>::subEntityTypes_{
      {RefElType::QUAD},{RefElType::SEGMENT,RefElType::SEGMENT,RefElType::SEGMENT,RefElType::SEGMENT}, {RefElType::POINT,RefElType::POINT,RefElType::POINT,RefElType::POINT}};

#if !defined (__INTEL_COMPILER)
    template<>
#endif
    const std::vector<iiVec_t> ReferenceElement<RefElType::QUAD>::subEntityCorners_{ 
      {{0,1,2,3}}, {{0,1}, {1,2}, {2,3}, {3,0}},{{0},{1},{2},{3}}};

#if !defined (__INTEL_COMPILER)
    template<>
#endif
    const fixedMatrix_t<2,4> ReferenceElement<RefElType::QUAD>::cornerCoord_(
                                                                             (fixedMatrix_t<2,4>() << 0,1,1,0,
                                                                              0,0,1,1).finished());

    //--------------------------------------------------------------------------
    // TETRA
    //--------------------------------------------------------------------------
#if !defined (__INTEL_COMPILER)
    template<>
#endif
    const itVec_t ReferenceElement<RefElType::TETRA>::subEntityTypes_ {
      {RefElType::TETRA},
        {RefElType::TRIA,RefElType::TRIA,RefElType::TRIA,RefElType::TRIA},
          {RefElType::SEGMENT,RefElType::SEGMENT,RefElType::SEGMENT,RefElType::SEGMENT,
              RefElType::SEGMENT,RefElType::SEGMENT},
            {RefElType::POINT,RefElType::POINT,RefElType::POINT,RefElType::POINT}
    };

#if !defined (__INTEL_COMPILER)
    template<>
#endif
    const std::vector<iiVec_t> ReferenceElement<RefElType::TETRA>::subEntityCorners_{
      {{0,1,2,3}},{{0,2,1}, {0,3,2}, {0,1,3}, {1,2,3}},
                    {{0,1}, {0,2}, {0,3}, {1,2}, {1,3}, {2,3}},{{0},{1},{2},{3}}};
#if !defined (__INTEL_COMPILER)
    template<>
#endif
    const fixedMatrix_t<3,4> ReferenceElement<RefElType::TETRA>::cornerCoord_
    ((fixedMatrix_t<3,4>() << 
      0,1,0,0,
      0,0,1,0,
      0,0,0,1).finished());

    //--------------------------------------------------------------------------
    // HEX
    //--------------------------------------------------------------------------
#if !defined (__INTEL_COMPILER)
    template<>
#endif
    const itVec_t ReferenceElement<RefElType::HEXA>::subEntityTypes_{
      {RefElType::HEXA},{RefElType::QUAD,RefElType::QUAD,RefElType::QUAD,RefElType::QUAD,RefElType::QUAD,RefElType::QUAD},{RefElType::SEGMENT,RefElType::SEGMENT,RefElType::SEGMENT,RefElType::SEGMENT,
                            RefElType::SEGMENT,RefElType::SEGMENT,RefElType::SEGMENT,RefElType::SEGMENT,RefElType::SEGMENT,RefElType::SEGMENT,RefElType::SEGMENT,RefElType::SEGMENT},
                                                                                                                            {RefElType::POINT,RefElType::POINT,RefElType::POINT,RefElType::POINT,RefElType::POINT,RefElType::POINT,RefElType::POINT,RefElType::POINT}
    };

#if !defined (__INTEL_COMPILER)
    template<>
#endif
    const std::vector<iiVec_t> ReferenceElement<RefElType::HEXA>::subEntityCorners_{
      {{0,1,2,3,4,5,6,7}},{{0,3,2,1}, {4,5,6,7}, {0,4,7,3}, {1,2,6,5}, {2,3,7,6},
                                                                         {0,1,5,4}},{{0,1},{0,3},{0,4},{1,2},{1,5},{2,3},{2,6},{3,7},{4,5},{4,7},
                                                                                                                                             {5,6},{6,7}},{{0},{1},{2},{3},{4},{5},{6},{7}}};
#if !defined (__INTEL_COMPILER)
    template<>
#endif
    const fixedMatrix_t<3,8> ReferenceElement<RefElType::HEXA>::cornerCoord_(
                                                                             (fixedMatrix_t<3,8>() << -1, 1, 1,-1,-1, 1, 1,-1,
                                                                              -1,-1, 1, 1,-1,-1, 1, 1,
                                                                              -1,-1,-1,-1, 1, 1, 1, 1).finished());

    //--------------------------------------------------------------------------
    // PRISM
    //--------------------------------------------------------------------------
#if !defined (__INTEL_COMPILER)
    template<>
#endif
    const itVec_t ReferenceElement<RefElType::PRISM>::subEntityTypes_{
      {RefElType::PRISM},{RefElType::TRIA,RefElType::TRIA,RefElType::QUAD,RefElType::QUAD,RefElType::QUAD},{RefElType::SEGMENT,RefElType::SEGMENT,RefElType::SEGMENT,RefElType::SEGMENT,
                             RefElType::SEGMENT,RefElType::SEGMENT,RefElType::SEGMENT,RefElType::SEGMENT,RefElType::SEGMENT},{RefElType::POINT,RefElType::POINT,RefElType::POINT,RefElType::POINT,RefElType::POINT,RefElType::POINT}};

#if !defined (__INTEL_COMPILER)
    template<>
#endif
    const std::vector<iiVec_t> ReferenceElement<RefElType::PRISM>::subEntityCorners_{
      {{0,1,2,3,4,5}},{{0,2,1}, {3,4,5}, {1,2,5,4}, {0,3,5,2}, {0,1,4,3}},
                        {{0,1},{0,2}, {0,3}, {1,2}, {1,4}, {2,5}, {3,4}, {3,5}, {4,5}},
                          {{0},{1},{2},{3},{4},{5}}};

#if !defined (__INTEL_COMPILER)
    template<>
#endif
    const fixedMatrix_t<3,6> ReferenceElement<RefElType::PRISM>::cornerCoord_
    ((fixedMatrix_t<3,6>() << 0 , 1, 0,0 ,1 ,0 ,
      0 , 0, 1,0 ,0 ,1 ,
      -1,-1,-1,1 ,1 ,1 ).finished());

    //--------------------------------------------------------------------------
    // PYRAMID
    //--------------------------------------------------------------------------
#if !defined (__INTEL_COMPILER)
    template<>
#endif
    const itVec_t ReferenceElement<RefElType::PYRAMID>::subEntityTypes_{
      {RefElType::PYRAMID},{RefElType::QUAD,RefElType::TRIA,RefElType::TRIA,RefElType::TRIA,RefElType::TRIA},{RefElType::SEGMENT,RefElType::SEGMENT,RefElType::SEGMENT,RefElType::SEGMENT,
                               RefElType::SEGMENT,RefElType::SEGMENT,RefElType::SEGMENT,RefElType::SEGMENT},{RefElType::POINT,RefElType::POINT,RefElType::POINT,RefElType::POINT,RefElType::POINT}};

#if !defined (__INTEL_COMPILER)
    template<>
#endif
    const std::vector<iiVec_t> ReferenceElement<RefElType::PYRAMID>::subEntityCorners_{
      {{0,1,2,3,4}},{{0,3,2,1}, {0,1,4}, {1,2,4}, {2,3,4}, {0,4,3}},{{0,1}, {0,3},
                                                                              {0,4}, {1,2}, {1,4}, {2,3}, {2,4}, {3,4}},{{0},{1},{2},{3},{4}}};

#if !defined (__INTEL_COMPILER)
    template<>
#endif
    const fixedMatrix_t<3,5> ReferenceElement<RefElType::PYRAMID>::cornerCoord_
    ((fixedMatrix_t<3,5>() << -1, 1, 1,-1, 0,
      -1,-1, 1, 1, 0,
      0, 0, 0, 0, 1).finished());
  } // end namespace base
} // end namespace eth

#endif // ETH_BASE_REF_EL_TYPES_I_HPP
