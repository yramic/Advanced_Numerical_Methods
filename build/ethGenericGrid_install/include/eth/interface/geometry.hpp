/********************************************************************
  created:	2013/06/27
  created:	27:6:2013   16:00
  filename: 	geometry.hpp
  author:		Raffael Casagrande
  
  purpose:	Provide abstract interface for a Geometry object.
*********************************************************************/

#ifndef _HPP_ACC43259_552B_4C80_8427_F88682275B54
#define _HPP_ACC43259_552B_4C80_8427_F88682275B54

namespace eth{
  namespace grid{
    template<class GRID_TRAITS, class GEOM_IMPL>
    class Geometry;
  }
}

// Implementation
//////////////////////////////////////////////////////////////////////////

#include "eth_base/CRTP_check.hpp"
#include "storage_types.hpp"
#include "eth_base/ref_el_types.hpp"

/**
 * \brief	Provides the mapping of an Entity/Intersection from a reference
 * 			element to the actual element in physical coordinates or some other
 *      reference element.
 * \tparam GRID_TRAITS a class which models the \ref GridTraitsDoc
 * \tparam GEOM_IMPL the underlying implementation class which is wrapped by this
 * 					 interface (engine pattern).
 * 					 
 * A Geometry object defines a map \f$ \Phi : R \mapsto W \f$ where 
 * \f$ R \subseteq \mathbb{R}^\text{dimFrom} \f$ and 
 * 	   \f$ W \subseteq \mathbb{R}^\text{dimTo} \f$. Here \f$ R \f$ refers to 
 * a generic reference element which is one of the list specified in 
 * eth::base::RefElType enum.
 * 
 * \note A Geometry Object is neither copyable nor assignable. The only way it
 * 		 can be passed is by using the Move constructor (or passing a
 * 		 reference). As soon as the GeometryObject is destroyed, the underlying
 * 		 implementation may also be destroyed and memory is released
 * 		 (depends on storage pattern, see below).
 * \note A Geometry object is used by many other classes. Entities use it to define
 *       a mapping from the local reference element to the physical element 
 *       (this is specified through gridTraits_t::geometry_t<codim>) as well as
 *       a mapping which describes how a child element is geometrically embedded
 *       in its father-element (gridTraits_t::localGeometry_t<dimFrom,dimTo>).
 *       Intersections use it to define a mapping from the intersection's reference
 *       element to the adjacent elements' reference elements
 *       (gridTraits_t::localGeometry_t<dimFrom,dimTo>). 
 *       The template parameter GEOM_IMPL of this class determines the type of
 *       this mapping.
 *
 * 			
 * <H3> Information for Grid-implementor </H3>
 * The Geometry implementation must provide two constants:
 *
 *   static constant          | explanation
 *   -------------------------|------------------------------------------------
 *   dimFrom                  | dimension of the space from which the geometry object maps.
 *   dimTo                    | dimension of the space to which this geometry object maps.
 *
 * - The Geometry Wrapper supports two types for storing the internal  
 *   implementation. The appropriate technique is determined based on the
 *   constant storageType which must be implemented in the engine:
 *   - if storageType = StorageType::UniquePtr, the wrapper uses internally  
 *	   a std::unique_ptr. In order to construct a wrapper, the implementor must
 *	   pass a std::unique_ptr<GEOM_IMPL> to the constructor of this wrapper. 
 *	   As soon as the user destructs the GeometryObject the underyling
 *	   implementation object (pointed to by unique_ptr) is also destroyed and
 *	   memory is released. This gives the user the flexibility to choose when
 *	   memory should be released.
 *	 - if storageType=StorageType::ConstRef the wrapper class stores a const  
 *	   REFERENCE to the underlying implementation class. This implies that the
 *	   reference passed in the constructor of this class  
 *     must remain valid for at least as long as the related Entity/Intersection
 *     remains valid. (In this scenario the Geometry implementation should be
 *     stored as a part of the entity implementation.) See also \ref LifetimeDoc
 *     .
 */
template<class GRID_TRAITS, class GEOM_IMPL>
class eth::grid::Geometry {
public:
  /// A model of GridTraits (see \ref GridTraitsDoc)
  typedef GRID_TRAITS gridTraits_t;
  /// Underlying implementation class (Engine Pattern)
  typedef GEOM_IMPL impl_t;

  typedef typename gridTraits_t::size_type size_type;

  /// dimension from which we map points.
  static const int dimFrom = impl_t::dimFrom;
  /// dimension to which the points are mapped.
  static const int dimTo   = impl_t::dimTo;

  /// type of global coordinates.
  typedef typename gridTraits_t::template fixedSizeMatrix_t<dimTo,1> globalCoord_t;
  /// type of reference element local coordinates.
  typedef typename gridTraits_t::template fixedSizeMatrix_t<dimFrom,1> localCoord_t;

  static const StorageType storageType = impl_t::storageType;
  // TODO: Uncomment this line as soon as it is supported by Intel Compiler.
  /*static_assert(storageType == StorageType::UniquePtr || 
    storageType == StorageType::ConstRef, 
    "Only Storage Types UniquePtr and ConstRef are supported by Geometry");*/


  //////////////////////////////////////////////////////////////////////////

  /// Forbid assignment
  Geometry& operator=(const Geometry& rhs) = delete;

  // Forbid copy construction:
  Geometry(const Geometry& rhs) = delete;

  /// Implementation Constructor.
  Geometry(typename StorageHelper<storageType,impl_t>::storageType_t impl) : 
    impl_(std::move(impl)) {}

  /// Move constructor:
  Geometry(Geometry&& other) : impl_(std::move(other.impl_)) {}

  /// Is the transformation described by this Geometry object affine?
  bool isAffine() const {
    return StorageHelper<storageType,impl_t>::ref(impl_).isAffine();
  }

  /// Type of reference element?
  eth::base::RefElType refElType() const {
    ETH_ASSERT_MSG(eth::base::ReferenceElements::getDimension(
      StorageHelper<storageType,impl_t>::ref(impl_).refElType())==dimFrom,
      "dim of Reference element != dimFrom");
    return StorageHelper<storageType,impl_t>::
      ref(impl_).refElType();
  }

  /**
   * \brief	Gets the number corners of reference element.
   * \remark  This method is exactly the same as 
   * 			ReferenceElement::numCorners(), it is here for convenience.
   */
  size_type numCorners() const {
    return StorageHelper<storageType,impl_t>::
      ref(impl_).numCorners();
  }

  /**
   * \brief	Gets the total number of nodes this geometry object contains
   */
  size_type numNodes() const {
    return StorageHelper<storageType,impl_t>::
      ref(impl_).numNodes();
  }


  /// Global coordinate of corner i
  globalCoord_t mapCorner(size_type i) const {
    return StorageHelper<storageType,impl_t>::
      ref(impl_).mapCorner(i);
  }

  /**
  * \brief	Map a list of points expressed in local coordinates to global
  * 			coordinates. (coordinates are stored as column vectors)
  * \param local A matrix of size dimFrom x NUM_POINTS.
  * \tparam NUM_POINTS the number of points to be mapped.
  */
  template<int NUM_POINTS>
  typename gridTraits_t::template fixedSizeMatrix_t<dimTo,NUM_POINTS> global
  (const typename gridTraits_t::template fixedSizeMatrix_t<dimFrom,NUM_POINTS>& local) const {
      return StorageHelper<storageType,impl_t>::
        ref(impl_).global(local);
  }

  /**
  * \brief	Map a list of points in global coordinates to the local 
  * 			coordinates of the reference element. (coordinates are stored
  * 			as column vectors).
  * \param global A matrix of size dimTo x NUM_POINTS
  * \tparam NUM_POINTS the number of point to be mapped.
  */
  template<int NUM_POINTS>
  typename gridTraits_t::template fixedSizeMatrix_t<dimFrom,NUM_POINTS> local
    (const typename gridTraits_t::template fixedSizeMatrix_t<dimTo,NUM_POINTS>& global) const {
      return StorageHelper<storageType,impl_t>::
        ref(impl_).local(global);
  }

  /**
   * \brief	Return the integration element at a set of points (in local
   * 			coordinates) (factor appearing in integral transformation
   * 			formula)
   * \tparam	NUM_POINTS the number of points at which the integration element
   * 					   should be evaluated.
   * \param	local	A matrix of size dimFromxNUM_POINTS
   * 					
   * Let \f$ \Phi : R \mapsto W \f$ denote the transformation represented by
   * this Geometry object and let \f$ J_\Phi(\xi) \f$ denote its Jacobian
   * (see Geometry::jacobianTransposed() ). The integration element \f$ g\f$
   * at point \f$ \xi \in R \f$ (local coordinates) is then given by
   * \f[
   *    g(\xi) = \sqrt{\left|\operatorname{det} J_\Phi^T(\xi) J_\Phi(\xi)\right|}.
   * \f]
   */
  template<int NUM_POINTS>
  typename gridTraits_t::template fixedSizeMatrix_t<1,NUM_POINTS> 
  integrationElement
  (const typename gridTraits_t::template fixedSizeMatrix_t<dimFrom,NUM_POINTS>& local) 
    const 
  {
      return StorageHelper<storageType,impl_t>::
        ref(impl_).integrationElement(local);
  }

/**
* \brief Get the integration element at the quadrature points defined by the
*        given quadrature rule...
* \tparam QUAD_RULE A class which fulfills the \ref QuadRuleDoc .
* \note  This method may have better performance than the 
*        integrationElement(const typename gridTraits_t::template fixedSizeMatrix_t< dimFrom, NUM_POINTS > &local) const
*        method
*        because the quadrature points are already available when the application
*        is launched (they are static). Thus the values of the shape functions
*        at the quadrature points must be calculated only once.
*/
  template<class QUAD_RULE>
  typename gridTraits_t::template fixedSizeMatrix_t<1,QUAD_RULE::numPoints>
     integrationElement()  const {
       static_assert(dimFrom == eth::base::ReferenceElement<QUAD_RULE::refElType>::dimension,
            "dimension of Reference Element doesn't match dimension of QuadRule");
      return StorageHelper<storageType,impl_t>::ref(impl_).template
        integrationElement<QUAD_RULE>();
  }

  /// The volume of the mapped element (not of the reference element!)
  typename gridTraits_t::ctype_t volume() const {
    return StorageHelper<storageType,impl_t>::ref(impl_).volume();
  }

  
  /**
   * \brief	Return the center of the mapped element.
   * \note	Depending on the implementation this can either be the centroid
   * 			or the centroid of all the mapped points or the mapped centroid
   * 			of the reference element etc..
   */
  globalCoord_t center() const {
    return StorageHelper<storageType,impl_t>::ref(impl_).center();
  }

  /**
  * \brief	Return one or more transposed jacobian(s) of the mapping at
  * 			\f$ \Phi : R \mapsto W \f$ at a given point(s)
  * 			\f$ \xi \in R \f$.
  * \tparam	NUM_POINTS	The number of points at which the jacobian should
  * 						be evaluated.
  * \param	local	The point(s) at which the jacobian should be returned.
  * 					Matrix of size dimFrom x NUM_POINTS, each column vector
  * 					corresponds to a point.
  * \returns A list of \f$ J_\Phi^T \f$. If NUM_POINTS>1, the returned matrix
  * 		   contains the \e transposed jacobian matrices as columns.
  * 
  * Let \f$ \Phi : R \mapsto W \f$ denote a mapping that maps the Reference
  * element D to its global representation (contained in this geometry object).
  * The Jacobian \f$ J_\Phi \f$ is the given by the  \f$ \text{dimTo } \times
  * \text{dimFrom} \f$
  * matrix
  * \f[
  *    J_\Phi(\xi) = \begin{pmatrix}
  *    \frac{\partial \Phi_1}{\partial \xi_1} & \cdots & 
  *    				\frac{\partial \Phi_1}{\partial \xi_{\text{dimFrom}}} \\
  *    \vdots & \ddots & \vdots \\
  *    \frac{\partial \Phi_{\text{dimTo}}}{\partial \xi_1} & \cdots & 
  *    \frac{\partial \Phi_{\text{dimTo}}}{\partial \xi_{\text{dimFrom}}}
  *    \end{pmatrix}
  * \f]
  */
  template<int NUM_POINTS>
  typename gridTraits_t::template fixedSizeMatrix_t<dimFrom,
    dimTo*NUM_POINTS> 
    jacobianTransposed(const typename gridTraits_t::template 
    fixedSizeMatrix_t<dimFrom,NUM_POINTS>& local) const {
      return StorageHelper<storageType,impl_t>::
        ref(impl_).jacobianTransposed(local);
  }

  /**
  * \brief Return the Jacobian at the quadrature points defined by QUAD_RULE.
  * \tparam QUAD_RULE   A class which fulfills the \ref QuadRuleDoc . 
  * \return A Matrix of size dimFrom x (QUAD_RULE::numPoints*dimTo) which
  *         contains the transposed jacobians (see 
  *         jacobianTransposed(const typename gridTraits_t::template fixedSizeMatrix_t< dimFrom, NUM_POINTS > &local) const
  *         ) as "column matrices" at the quadrature points defined by QUAD_RULE.
  */
  template<class QUAD_RULE>
  typename gridTraits_t::template fixedSizeMatrix_t<dimFrom,
      QUAD_RULE::numPoints*dimTo> jacobianTransposed() const {
       static_assert(dimFrom == eth::base::ReferenceElement<QUAD_RULE::refElType>::dimension,
            "dimension of Reference Element doesn't match dimension of QuadRule");
          return StorageHelper<storageType,impl_t>::ref(impl_).template
              jacobianTransposed<QUAD_RULE>();
  }

  /**
  * \brief	Return the inverse of the transposed Jacobian,
  * 			\f$ J_\Phi^{-T} \f$. See jacobianTransposed() for more details.
  * \tparam	NUM_POINTS the number of points at which the Jacobian is evaluated.
  * \note		if dimFrom < dimTo, then this method returns the inverse of
  *         the gram Matrix \f$ G := J_\Phi(\xi)^T J_\Phi(\xi) \f$.
  */
  template<int NUM_POINTS>
  typename gridTraits_t::template fixedSizeMatrix_t<dimFrom,
    dimFrom*NUM_POINTS> 
    jacobianInverseTransposed(const typename gridTraits_t::template 
    fixedSizeMatrix_t<dimFrom,NUM_POINTS>& local) const {
      return StorageHelper<storageType,impl_t>::
        ref(impl_).jacobianInverseTransposed(local);
  }

  /**
  * \brief Return the inverse of the transposed Jacobians at the quadrature
  *        points defined by QUAD_RULE
  * \tparam QUAD_RULE A class which fulfills the \ref QuadRuleDoc concept.
  * \note		if dimFrom < dimTo, then this method returns the inverse of
  *         the gram Matrix \f$ G := J_\Phi(\xi)^T J_\Phi(\xi) \f$.
  */
  template<class QUAD_RULE>
  typename gridTraits_t::template fixedSizeMatrix_t<dimFrom,
      QUAD_RULE::numPoints*dimFrom> jacobianInverseTransposed() const {
       static_assert(dimFrom == eth::base::ReferenceElement<QUAD_RULE::refElType>::dimension,
            "dimension of Reference Element doesn't match dimension of QuadRule");
      return StorageHelper<storageType,impl_t>::ref(impl_).template
          jacobianInverseTransposed<QUAD_RULE>();
  }
    
  /**
   * \brief Return a constant reference to the underlying geometry implementation
   * \note This method should be used only(!) when the interface methods do not fit!
   */
  const impl_t& impl( ) const { return StorageHelper<storageType,impl_t>::ref(impl_); }

private:
  /// Underlying implementation...
  typename StorageHelper<storageType,impl_t>::storageType_t	 impl_;
};




#endif // _HPP_ACC43259_552B_4C80_8427_F88682275B54
