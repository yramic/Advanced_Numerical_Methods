/***********************************************************************
 *                                                                     *
 * Code for Course "Advanced Numerical Methods for CSE"                *
 * (Prof. Dr. R. Hiptmair)                                             * 
 * Author:                                                             *
 * Date:                                                               *
 * (C) Seminar for Applied Mathematics, ETH Zurich                     *
 * This code can be freely used for non-commercial purposes as long    *
 * as this header is left intact.                                      *
 ***********************************************************************/
#ifndef KERNEL_HPP
#define KERNEL_HPP


/**
* \brief Kernel functor $\frac{num}{|x-y|}$ if $x != y$, else 0
*/
class Kernel
{
public:

    /**
    * \brief Constructor
    */
    Kernel(double num):
        num_(num)
    { }

    /**
    * \brief Functor
    */
    double operator()(double x, double y);

private:

    double num_; // numerator
};

#endif // KERNEL_HPP
