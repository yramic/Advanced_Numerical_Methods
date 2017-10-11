#ifndef IS_ADMISSIBLE_HPP
#define IS_ADMISSIBLE_HPP

#include "../include/node.hpp"

/**
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
    double get_max();
    double get_min();
    bool is_admissible();
};

class AdmissibilityD: public Admissibility
{
public:
    double get_max(double xl, double xr, double yl, double yr);
    double get_min(double xl, double xr, double yl, double yr);
    bool is_admissible(double xl, double xr, double yl, double yr, double eta);

};

class AdmissibilityH: public Admissibility
{
public:
    double get_max(Node* x, Node* y);
    double get_min(Node* x, Node* y);
    bool is_admissible(Node* x, Node* y, double eta);

};



#endif // IS_ADMISSIBLE_HPP
