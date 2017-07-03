#ifndef KERNEL_HPP
#define KERNEL_HPP


class Kernel
{
public:

    Kernel( double num ):
        num_( num )
    { }

    double operator()( double x, double y );

private:

    double num_;
};

#endif // KERNEL_HPP
