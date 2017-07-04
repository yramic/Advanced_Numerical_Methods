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
