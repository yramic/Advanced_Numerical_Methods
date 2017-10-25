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

    /*!
    * \brief Virtual Functor
    * \param x1 x coordinate of first point
    * \param y1 y coordinate of first point
    * \param x2 x coordinate of second point
    * \param y2 y coordinate of second point
    */
    virtual double operator()(double x1, double y1, double x2, double y2) = 0;

protected:
    double num_; //!< numerator
};


/*!
* \brief Kernel functor \f$-\frac{1}{2 \pi}\log{\left|\vec{x}-\vec{y}\right|}\f$
*/
class KernelGalerkin: public Kernel
{
public:

    /*!
    * \brief Functor for 2D
    * \param x1 x coordinate of first point
    * \param y1 y coordinate of first point
    * \param x2 x coordinate of second point
    * \param y2 y coordinate of second point
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
    * \brief Functor for 2D Singular Kernel
    * \param x1 x coordinate of first point
    * \param y1 y coordinate of first point
    * \param x2 x coordinate of second point
    * \param y2 y coordinate of second point
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
    * \brief Functor for 2D Singular Kernel
    * \param x1 x coordinate of first point
    * \param y1 y coordinate of first point
    * \param x2 x coordinate of second point
    * \param y2 y coordinate of second point
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
    * \brief Functor for 2D
    * \param x1 x coordinate of first point
    * \param y1 y coordinate of first point
    * \param x2 x coordinate of second point
    * \param y2 y coordinate of second point
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
    * \brief Functor for 2D
    * \param x1 x coordinate of first point
    * \param y1 y coordinate of first point
    * \param x2 x coordinate of second point
    * \param y2 y coordinate of second point
    */
    double operator()(double x1, double y1, double x2, double y2);
};

/**
* \brief Kernel functor \f$\cos(|\vec{x}-\vec{y}|)\f$
*/
class GlobalSmoothKernel: public Kernel
{
public:
    /**
    * \brief Functor for 2D
    * \param x1 x coordinate of first point
    * \param y1 y coordinate of first point
    * \param x2 x coordinate of second point
    * \param y2 y coordinate of second point
    */
    double operator()(double x1, double y1, double x2, double y2);
};

/**
* \brief Kernel functor \f$e^{-|\vec{x}-\vec{y}|^2}\f$
*/
class GaussKernel: public Kernel
{
public:
    /**
    * \brief Functor for 2D
    * \param x1 x coordinate of first point
    * \param y1 y coordinate of first point
    * \param x2 x coordinate of second point
    * \param y2 y coordinate of second point
    */
    double operator()(double x1, double y1, double x2, double y2);
};
#endif // KERNEL_HPP

