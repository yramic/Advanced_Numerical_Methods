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
    Kernel() {}
    Kernel(double num):
        num_(num)
    { }

    /**
    * \brief Functor
    */
    //double operator()(double x, double y);
    //double operator()(double x1, double y1, double x2, double y2);
    double operator()();
protected:

    double num_; // numerator
};

class Kernel2D: public Kernel
{
public:

    /**
    * \brief Constructors
    */
    Kernel2D(): Kernel() {}
    Kernel2D(double num): Kernel(num) {}

    /**
    * \brief Functor
    */
    double operator()(double x, double y);
};

class Kernel4D: public Kernel
{
public:

    /**
    * \brief Functor
    */
    double operator()(double x1, double y1, double x2, double y2);
};

class PolynomialKernel: public Kernel
{
public:

    /**
    * \brief Functor
    */
    double operator()(double x1, double y1, double x2, double y2);
};

#endif // KERNEL_HPP
