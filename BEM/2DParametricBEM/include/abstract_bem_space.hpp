/**
 * \file abstract_bem_space.hpp
 * \brief This file declares an abstract class representing a BEM Space
 *
 * This File is a part of the 2D-Parametric BEM package
 */

#ifndef ABSTRACTBEMSPACEHPP
#define ABSTRACTBEMSPACEHPP

#include <vector>
#include <utility>

namespace parametricbem2d {
  /**
   * \class AbstractBEMSpace
   * \brief This abstract class declares the interface for a BEM Space
   */
  class AbstractBEMSpace {
  public:
    /**
     * This typedef aids in making a std::vector of function pointers.
     * The basis functions have the signature: double basis_function(double);
     */
    typedef double(*BasisFunctionPointer)(double t);

    /**
     * This typedef aids in returning a function pointer for local to global map
     * defined in 1.4.75. The maps have the following signature:
     * unsigned int map (unsigned int q, unsigned int n, unsigned int N);
     *
     * @param q Index of the local/reference shape function
     * @param n Index of the panel for which this map is applied
     * @param N Total number of panels in the system
     * @return Integer denoting the global shape function number corresponding
     *         to the local shape function referenced by q for the panel no. n
     */
    typedef unsigned int(*LocGlobMapPointer) (unsigned int q, unsigned int n, unsigned int N);

    /**
      * This function is used for querying the parameter interval.
      * The standard parameter interval [-1,1] is used and it can't
      * be overriden in the inherited classes as the function is non
      * virtual. The function is declared static as it is independent
      * of the concrete object.
      *
      * @return A std::pair<double,double> object containing
      *          the valid parameter range for parametrization
      *          that is [-1,1]
      */
     static std::pair<double,double> ParameterRange(void) {
       // Parameter range : [-1,1]
       return std::make_pair(-1.,1.);
     }

    /**
     * This function is used to get the value of Q for a BEM Space.
     * Q stands for the number of local shape functions for the Space.
     *
     * @return Integer Q which is the number of local shape functions
     */
    int getQ() const {
      return q_;
    }

    /**
     * This function is used to get a vector containing all the Reference
     * Shape Functions for the given BEM Space.
     *
     * @return An std::vector object of size Q containing all the Reference
     *         Shape Functions in the form of function pointers
     */
     std::vector<BasisFunctionPointer> getShapeFunctions() const {
       return referenceshapefunctions_;
     }

     /**
     * This function is used for checking whether a value t is within the
     * valid parameter range. This function is non virtual to prevent it
     * from being overriden as the parameter interval is fixed. It is
     * declared static because the check is independent of the concrete
     * implementation.
     *
     * @param t The value to be checked
     * @return boolean indicating result of the performed check
     */
    static bool IsWithinParameterRange(double t) {
        double a,b;
        std::tie(a,b) = ParameterRange();
        return ( t>=a && t<=b );
    }

    /**
    * This function returns a function pointer to the
    * locap to global map function
    *
    * @return function pointer to local to global map
    */
    LocGlobMapPointer getLocGlobMap() const{
      return locglobmap_;
    }

    /**
    * This function returns the dimension of the given BEM space on the given
    * number of panels. It is a pure virtual function and has to be implemented
    * by the inherited class.
    *
    * @param numpanels Number of panels in the mesh
    * @return Dimension of the BEM space
    */
    virtual unsigned int getSpaceDim(unsigned int numpanels) const=0;

   protected:
     /**
      * protected constructor ensures this base class is non instantiable.
      */
     AbstractBEMSpace():q_(0) {};
     /**
      * A protected vector containing the Reference Shape Functions.
      */
     std::vector<BasisFunctionPointer> referenceshapefunctions_;
     /**
      * A protected pointer of the type LocGlobMapPointer which is returned by
      * the function getLocGlobMap
      */
     LocGlobMapPointer locglobmap_;
     /**
      * A protected integer q which stores the number of reference shape
      * functions in derived concrete classes
      */
     int q_;
  }; // class AbstractBEMSpace
} // namespace parametricbem2d

#endif //ABSTRACTBEMSPACEHPP
