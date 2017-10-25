#ifndef IS_ADMISSIBLE_HPP
#define IS_ADMISSIBLE_HPP

#include "../include/node.hpp"

/*!
* \brief Primitive functions used to check whether a cluster is admissible
* (eta-admissibility)
*/

//double get_max(double xl, double xr, double yl, double yr);
//double get_min(double xl, double xr, double yl, double yr);
//bool is_admissible(double xl, double xr, double yl, double yr, double eta);


class Admissibility
{
public:
    /*
    virtual double get_max(double xl, double xr, double yl, double yr);
    virtual double get_min(double xl, double xr, double yl, double yr);
    virtual bool is_admissible(double xl, double xr, double yl, double yr, double eta);
    virtual double get_max(Node* a, Node* b);
    virtual double get_min(Node* a, Node* b);
    virtual bool is_admissible(Node* x, Node* y, double eta);
    */  
    /*!
     * \brief return maximum base
     */
    double get_max();
    /*!
     * \brief return minimum base
     */
    double get_min();
    /*!
     * \brief return if the cluster is admissible
     */
    bool is_admissible();
};
/*!
* \brief Class for 2D admissibility problems
*/
class AdmissibilityD: public Admissibility
{
public:
    /*!
     * \brief returns biggest edge
     */
    double get_max(double xl, double xr, double yl, double yr);
    /*!
     * \brief returns smallest edge
     */
    double get_min(double xl, double xr, double yl, double yr);
    /*!
     * \brief return if the cluster is admissible
     */
    bool is_admissible(double xl, double xr, double yl, double yr, double eta);

};
/*!
* \brief Class for 4D admissibility problems(admissibility between 2 bounding boxes)
*/
class AdmissibilityH: public Admissibility
{
public:
    /*!
     * \brief returns the biggest edge of the two bounding boxes
     */
    double get_max(Node* x, Node* y);
    /*!
     * \brief returns the distance of the two bounding boxes
     */
    double get_min(Node* x, Node* y);
    /*!
     * \brief returns if the two bounding boxes are admissible
     */
    bool is_admissible(Node* x, Node* y, double eta);

};



#endif // IS_ADMISSIBLE_HPP
