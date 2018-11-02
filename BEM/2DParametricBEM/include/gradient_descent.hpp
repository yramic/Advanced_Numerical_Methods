#ifndef GRADIENTDESCENT
#define GRADIENTDESCENT

#include <Eigen/Dense>
#include <cmath>
#include <cstdlib>

//srand(1);

template <typename T>
class helper{
public:
  static double distance(T x0, T x1) {
    return (x0-x1).norm();
  }
  static double minimum(T x0) {
    return x0.minCoeff();
  }
  static double maximum(T x0) {
    return x0.maxCoeff();
  }
  static Eigen::VectorXd RandomStart(T x0) {
    unsigned int N = x0.rows();
    return Eigen::VectorXd::Random(N);
  }
};

template<>
class helper<double> {
public:
  static double distance(double x0, double x1) {
    return fabs(x0-x1);
  }
  static double minimum(double x0) {
    return x0;
  }
  static double maximum(double x0) {
    return x0;
  }
  static double RandomStart(double x0) {
    return 1-2*rand()/RAND_MAX;
  }
};

template <typename Grad, typename Solution_t>
std::pair<Solution_t,bool> GradientDescent(Grad gradient, Solution_t initial) {
  Solution_t x0;
  Solution_t x1 = initial;
  unsigned int max_count = 10000;
  unsigned int count = 0;
  unsigned int max_random_starts = 5;
  unsigned int random_starts = 0;
  bool exists = true;
  do {
    ++count;
    x0 = x1;
    x1 = x1 - 0.1 * gradient(x1);
    // Point goes beyond boundary, random start
    if (helper<Solution_t>::minimum(x1) < -1 || helper<Solution_t>::maximum(x1) > 1) {
      std::cout << "random start as solution reached : " << x1 << std::endl;
      count = 0;
      ++random_starts;
      x1 = helper<Solution_t>::RandomStart(x1);
      if (random_starts == max_random_starts)
        exists = false;
    }
  } while (helper<Solution_t>::distance(x0,x1)>1e-8 && exists==true && count < max_count);
  std::cout << "final count = " << count << std::endl;
  if (count == max_count)
    exists = false;
  return std::make_pair(x1,exists);
}



#endif // GRADIENTDESCENT
