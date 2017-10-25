#ifndef KERNEL_HPP
#define KERNEL_HPP

/*!
* \brief Kernel base class
*/
class Kernel
{
public:

    /*!
    * \brief Default Constructor
    */
    Kernel() {}
    /*!
    * \brief Default Constructor
    */
    Kernel(double num):
        num_(num)
    { }

    // functor for 2D
    double operator()(double x, double y){}

    /*!
    * \brief Virtual Functor
    */
    virtual double operator()(double x1, double y1, double x2, double y2) = 0;
    //double operator()();
protected:

    double num_; //!< numerator
};

/*!
* \brief Kernel functor \f$\frac{num}{|x-y|}\f$ if \f$x != y\f$, else 0
*/
class Kernel2D: public Kernel
{
public:

    /*!
    * \brief Constructor
    */
    Kernel2D(): Kernel() {}
    /*!
    * \brief Constructor for num
    */
    Kernel2D(double num): Kernel(num) {}

    /*!
    * \brief Functor for 2D
    */
    double operator()(double x, double y);
    /*!
    * \brief Definition for virtual functor
    */
    double operator()(double x1, double y1, double x2, double y2);

};

/*!
* \brief Kernel functor \f$-\frac{1}{2 \pi}\log{\left|\vec{x}-\vec{y}\right|}\f$
*/
class Kernel4D: public Kernel
{
public:

    /*!
    * \brief Functor for 4D
    */
    double operator()(double x1, double y1, double x2, double y2);
};

/*!
* \brief Kernel functor \f$\log{\left|\vec{x}-\vec{y}\right|}\f$
*/
class SingularKernel: public Kernel
{
public:

    /*!
    * \brief Functor for 4D Singular Kernel
    */
    double operator()(double x1, double y1, double x2, double y2);
};

/*!
* \brief Kernel functor \f$\frac{1}{\log{\left|\vec{x}-\vec{y}\right|}}\f$
*/
class SingularKernelf: public Kernel
{
public:

    /*!
    * \brief Functor for 4D Singular Kernel
    */
    double operator()(double x1, double y1, double x2, double y2);
};

/*!
* \brief Kernel functor \f${x1}\times{x2}\times{y1}\times{y2}\f$
*/
class PolynomialKernel: public Kernel
{
public:

    /*!
    * \brief Functor for 4D
    */
    double operator()(double x1, double y1, double x2, double y2);
};

/*!
* \brief Kernel functor 1
*/
class ConstantKernel: public Kernel
{
public:
    /*!
    * \brief Default Constructor
    */
    ConstantKernel(): Kernel() {}
    /*!
    * \brief Constructor for num
    */
    ConstantKernel(double num):
        Kernel(num)
    { }
    /*!
    * \brief Functor for 4D
    */
    double operator()(double x1, double y1, double x2, double y2);
    /*!
    * \brief Functor for 2D
    */
    double operator()(double x1, double y2);
};

/**
* \brief Kernel functor \f$\cos(|\vec{x}-\vec{y}|)\f$
*/
class GlobalSmoothKernel: public Kernel
{
public:
    /**
    * \brief Functor for 4D
    */
    double operator()(double x1, double y1, double x2, double y2);
    /*!
    * \brief Functor for 2D
    */
    double operator()(double x1, double y2);
};

/**
* \brief Kernel functor \f$e^{-|\vec{x}-\vec{y}|^2}\f$
*/
class GaussKernel: public Kernel
{
public:
    /**
    * \brief Functor for 4D
    */
    double operator()(double x1, double y1, double x2, double y2);
};
#endif // KERNEL_HPP

