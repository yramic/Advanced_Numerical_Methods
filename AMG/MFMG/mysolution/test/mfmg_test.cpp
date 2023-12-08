/**
 * @file mfmg_test.cc
 * @brief ADVNCSE homework MFMG test code
 * @author Bob Schreiner, Jörg Nick
 * @date November 2023
 * @copyright Developed at SAM, ETH Zurich
 */

#include "../mfmg.h"

#include <gtest/gtest.h>

namespace MFMG::test {
TEST(MFMG, applyGridOperator) {
  const unsigned l = 7;
  const unsigned M = std::pow(2, l) - 1;
  //c=1
  DirichletBVPMultiGridSolver solver(1, l);
  GridFunction u = Eigen::MatrixXd::Zero(M + 2, M + 2);
  GridFunction analytic_result = Eigen::MatrixXd::Zero(M + 2, M + 2);
  GridFunction result = Eigen::MatrixXd::Zero(M + 2, M + 2);
  const double h = 1.0 / (M + 1);
  for (unsigned int i = 0; i < M + 2; ++i) {
    for (unsigned int j = 0; j < M + 2; ++j) {
      // The function $u(x,y)=\sin(\pi x)\sin(\pi y)$ vanishes along
      // the boundary and is an eigenfunction of the Laplacian $\Delta$.
      u(i, j) = std::sin(i * h * M_PI) * std::sin(j * h * M_PI);
    }
  }
  analytic_result = (M_PI * M_PI * 2 + 1) * u;
  result = solver.applyGridOperator(l, u);
  ASSERT_NEAR((result / (h * h) - analytic_result).cwiseAbs().maxCoeff(), 0,
              1e-3);
}
TEST(MFMG, directSolve) {
  const unsigned l = 5;
  const unsigned M = std::pow(2, l) - 1;
  DirichletBVPMultiGridSolver solver(1.0, l);
  GridFunction phi = Eigen::MatrixXd::Zero(M + 2, M + 2);
  const double h = 1.0 / (M + 1);

  //Differentiable function that goes to zero on the boundary
  auto f = [](double x, double y) {
    return (M_PI * M_PI * 2 + 1.0) * std::sin(x * M_PI) * std::sin(y * M_PI);
  };

  for (unsigned int i = 1; i < M + 1; ++i) {
    for (unsigned int j = 1; j < M + 1; ++j) {
      phi(i, j) = f(i * h, j * h);
    }
  }
  GridFunction result = solver.directSolve(l, h * h * phi);
  ASSERT_NEAR((result - phi / (M_PI * M_PI * 2 + 1.0)).cwiseAbs().maxCoeff(), 0,
              1e-3);
}
TEST(MFMG, prolongate_and_restrict) {
  const unsigned l = 4;
  const unsigned M = std::pow(2, l) - 1;
  const double h = 1.0 / (M + 1);
  DirichletBVPMultiGridSolver solver(1.0, l);
  GridFunction phi = Eigen::MatrixXd::Zero(M + 2, M + 2);

  for (unsigned int i = 1; i < M + 1; ++i) {
    for (unsigned int j = 1; j < M + 1; ++j) {
      phi(i, j) = std::sin(i * h * M_PI) * std::sin(j * h * M_PI);
    }
  }
  phi.col(0) = Eigen::VectorXd::Zero(M + 2);
  phi.col(M + 1) = Eigen::VectorXd::Zero(M + 2);
  phi.row(0) = Eigen::VectorXd::Zero(M + 2);
  phi.row(M + 1) = Eigen::VectorXd::Zero(M + 2);
  GridFunction phi_prolong = solver.prolongate(l + 1, phi);
  GridFunction phi_restr = solver.restrict(l + 1, phi_prolong);
  double error = (phi - 0.25 * phi_restr).cwiseAbs().maxCoeff();
  ASSERT_NEAR(error, 0, 1e-2);
}

TEST(MFMG, sweepGaussSeidel) {
  const unsigned l = 2;
  const unsigned M = std::pow(2, l) - 1;
  DirichletBVPMultiGridSolver solver(0, l);
  GridFunction phi = Eigen::MatrixXd::Zero(M + 2, M + 2);
  const double h = 1.0 / (M + 1);

  //Differentiable function that goes to zero on the boundary
  auto f = [](double x, double y) {
    return std::sin(x * M_PI) * std::sin(y * M_PI);
  };

  for (unsigned int i = 1; i < M + 1; ++i) {
    for (unsigned int j = 1; j < M + 1; ++j) {
      phi(i, j) = f(i * h, j * h);
    }
  }
  GridFunction mu = Eigen::MatrixXd::Zero(M + 2, M + 2);
  unsigned iter = 0;
  unsigned max_iter = 30;
  unsigned L0 = 2;
  do {
    solver.sweepGaussSeidel(l, mu, phi);
    iter++;
  } while (solver.residual(l, mu, phi).cwiseAbs().maxCoeff() > 1e-9 &&
           iter < max_iter);
  ASSERT_NEAR(solver.residual(l, mu, phi).cwiseAbs().maxCoeff(), 0, 1e-8);
}


TEST(MFMG, twolevelgridIteration) {
  const unsigned l = 6;
  const unsigned M = std::pow(2, l) - 1;
  DirichletBVPMultiGridSolver solver(1, l);
  GridFunction phi = Eigen::MatrixXd::Zero(M + 2, M + 2);
  const double h = 1.0 / (M + 1);

  //Differentiable function that goes to zero on the boundary
  auto f = [](double x, double y) {
    return std::sin(x * M_PI) * std::sin(y * M_PI);
  };

  for (unsigned int i = 1; i < M + 1; ++i) {
    for (unsigned int j = 1; j < M + 1; ++j) {
      phi(i, j) = h * h * f(i * h, j * h);
    }
  }
  GridFunction mu = Eigen::MatrixXd::Zero(M + 2, M + 2);
  unsigned iter = 0;
  unsigned max_iter = 20;
  unsigned L0 = 3;
  do {
    mu = solver.multigridIteration(mu, phi, L0);
    iter++;
  } while (solver.residual(l, mu, phi).cwiseAbs().maxCoeff() > 1e-9 &&
           iter < max_iter);
  ASSERT_NEAR(solver.residual(l, mu, phi).cwiseAbs().maxCoeff(), 0, 1e-9);
  ASSERT_NEAR(
      (mu - phi / (h * h * (2 * M_PI * M_PI + 1))).cwiseAbs().maxCoeff(), 0,
      1e-5);
}
TEST(MFMG, twolevelgridIterationLshape) {
  const unsigned l = 6;
  const unsigned M = std::pow(2, l) - 1;
  DirichletBVPMultiGridSolver solver(1, l);
  GridFunction phi = Eigen::MatrixXd::Zero(M + 2, M + 2);
  const double h = 1.0 / (M + 1);

  //Differentiable function that goes to zero on the boundary
  auto f = [](double x, double y) {
    return std::sin(x * M_PI) * std::sin(y * M_PI);
  };

  for (unsigned int i = 1; i < M + 1; ++i) {
    for (unsigned int j = 1; j < M + 1; ++j) {
      phi(i, j) = h * h * f(i * h, j * h);
    }
  }
  GridFunction mu = Eigen::MatrixXd::Zero(M + 2, M + 2);
  unsigned iter = 0;
  unsigned max_iter = 20;
  unsigned L0 = 3;
  do {
    mu = solver.multigridIteration(mu, phi, L0);
    iter++;
  } while (solver.residual(l, mu, phi).cwiseAbs().maxCoeff() > 1e-9 &&
           iter < max_iter);
  ASSERT_NEAR(solver.residual(l, mu, phi).cwiseAbs().maxCoeff(), 0, 1e-9);
  ASSERT_NEAR(
      (mu - phi / (h * h * (2 * M_PI * M_PI + 1))).cwiseAbs().maxCoeff(), 0,
      1e-5);
}
TEST(MFMG, multigridIteration) {
  const unsigned l = 9;
  const unsigned M = std::pow(2, l) - 1;
  DirichletBVPMultiGridSolver solver(1, l);
  GridFunction phi = Eigen::MatrixXd::Zero(M + 2, M + 2);
  const double h = 1.0 / (M + 1);

  //Differentiable function that goes to zero on the boundary
  auto f = [](double x, double y) {
    return std::sin(x * M_PI) * std::sin(y * M_PI);
  };

  for (unsigned int i = 1; i < M + 1; ++i) {
    for (unsigned int j = 1; j < M + 1; ++j) {
      phi(i, j) = h * h * f(i * h, j * h);
    }
  }
  GridFunction mu = Eigen::MatrixXd::Zero(M + 2, M + 2);
  unsigned iter = 0;
  unsigned max_iter = 20;
  unsigned L0 = 2;
  do {
    mu = solver.multigridIteration(mu, phi, L0);
    iter++;
  } while (solver.residual(l, mu, phi).cwiseAbs().maxCoeff() > 1e-9 &&
           iter < max_iter);
  ASSERT_NEAR(solver.residual(l, mu, phi).cwiseAbs().maxCoeff(), 0, 1e-9);
  ASSERT_NEAR(
      (mu - phi / (h * h * (2 * M_PI * M_PI + 1))).cwiseAbs().maxCoeff(), 0,
      1e-5);
}
}  // namespace MFMG::test

namespace MFMGLshape::test {
TEST(MFMG, applyGridOperatorLshape) {
  const unsigned l = 12;
  const unsigned M = std::pow(2, l) - 1;
  //c=1
  DirichletBVPMultiGridSolver solver(1, l);
  GridFunction u = Eigen::MatrixXd::Zero(M + 2, M + 2);
  GridFunction analytic_result = Eigen::MatrixXd::Zero(M + 2, M + 2);
  GridFunction result = Eigen::MatrixXd::Zero(M + 2, M + 2);
  const double h = 1.0 / (M + 1);
  const double midpoint = (M + 1) / 2;
  for (unsigned int i = 0; i < M + 2; ++i) {
    for (unsigned int j = 0; j < M + 2; ++j) {
      // The function $u(x,y)=\sin(\pi x)\sin(\pi y)$ vanishes along
      // the boundary and is an eigenfunction of the Laplacian $\Delta$.
      u(i, j) = std::sin(i * h * 2 * M_PI) * std::sin(j * h * 2 * M_PI);
      if (i <= midpoint && j >= midpoint) {
        u(i, j) = 0;
      }
    }
  }
  result = solver.applyGridOperator(l, u);
  analytic_result = (M_PI * M_PI * 8 + 1) * u;
  ASSERT_NEAR((result / (h * h) - analytic_result).cwiseAbs().maxCoeff(), 0,
              1e-4);
}
TEST(MFMG, directSolveLshape) {
  const unsigned l = 3;
  const unsigned M = std::pow(2, l) - 1;
  DirichletBVPMultiGridSolver solver(1.0, l);
  GridFunction phi = Eigen::MatrixXd::Zero(M + 2, M + 2);
  const double h = 1.0 / (M + 1);

  //Differentiable function that goes to zero on the boundary
  auto f = [](double x, double y) {
    return std::sin(2 * x * M_PI) * std::sin(2 * y * M_PI);
  };
  int midpoint = (M + 1) / 2;
  for (unsigned int i = 1; i < M + 1; ++i) {
    for (unsigned int j = 1; j < M + 1; ++j) {
      phi(i, j) = f(i * h, j * h);
      if (i <= midpoint && j >= midpoint) {
        phi(i, j) = 0.0;
      }
    }
  }
  GridFunction result = solver.directSolve(l, h * h * phi);
  //std::cout<< "phi:"<<std::endl<< phi/(M_PI*M_PI*8+1.0)<<std::endl;
  //std::cout<< "result"<<std::endl<<result<<std::endl;
  ASSERT_NEAR((result - phi / (M_PI * M_PI * 8 + 1.0)).cwiseAbs().maxCoeff(), 0,
              1e-3);
}

TEST(MFMG, prolongate_and_restrict) {
  const unsigned l = 4;
  const unsigned M = std::pow(2, l) - 1;
  const double h = 1.0 / (M + 1);
  DirichletBVPMultiGridSolver solver(1.0, l);
  GridFunction phi = Eigen::MatrixXd::Zero(M + 2, M + 2);

  for (unsigned int i = 1; i < M + 1; ++i) {
    for (unsigned int j = 1; j < M + 1; ++j) {
      phi(i, j) = std::sin(i * h * M_PI) * std::sin(j * h * M_PI);
    }
  }
  phi.col(0) = Eigen::VectorXd::Zero(M + 2);
  phi.col(M + 1) = Eigen::VectorXd::Zero(M + 2);
  phi.row(0) = Eigen::VectorXd::Zero(M + 2);
  phi.row(M + 1) = Eigen::VectorXd::Zero(M + 2);
  GridFunction phi_prolong = solver.prolongate(l + 1, phi);
  GridFunction phi_restr = solver.restrict(l + 1, phi_prolong);
  double error = (phi - 0.25 * phi_restr).cwiseAbs().maxCoeff();
  ASSERT_NEAR(error, 0, 1e-2);
}

TEST(MFMG, sweepGaussSeidelLshape) {
  const unsigned l = 3;
  const unsigned M = std::pow(2, l) - 1;
  DirichletBVPMultiGridSolver solver(0, l);
  GridFunction phi = Eigen::MatrixXd::Zero(M + 2, M + 2);
  const double h = 1.0 / (M + 1);

  //Differentiable function that goes to zero on the boundary
  auto f = [](double x, double y) {
    return std::sin(x * 2 * M_PI) * std::sin(y * 2 * M_PI);
  };

  const double midpoint = (M + 1) / 2;
  for (unsigned int i = 1; i < M + 1; ++i) {
    for (unsigned int j = 1; j < M + 1; ++j) {
      phi(i, j) = f(i * h, j * h);
      if (i <= midpoint && j >= midpoint) {
        phi(i, j) = 0;
      }
    }
  }

  GridFunction mu = Eigen::MatrixXd::Zero(M + 2, M + 2);
  unsigned iter = 0;
  unsigned max_iter = 100;
  unsigned L0 = 4;
  do {
    solver.sweepGaussSeidel(l, mu, phi);
    iter++;
  } while (solver.residual(l, mu, phi).cwiseAbs().maxCoeff() > 1e-9 &&
           iter < max_iter);
  ASSERT_NEAR(solver.residual(l, mu, phi).cwiseAbs().maxCoeff(), 0, 1e-8);
}

TEST(MFMG, multilevelgridIterationLshape) {
  const unsigned l = 10;
  const unsigned M = std::pow(2, l) - 1;
  DirichletBVPMultiGridSolver solver(1, l);
  GridFunction phi = Eigen::MatrixXd::Zero(M + 2, M + 2);
  const double h = 1.0 / (M + 1);
  const double midpoint = (M + 1) / 2;
  //Differentiable function that goes to zero on the boundary
  auto f = [](double x, double y) {
    return std::sin(x * 2 * M_PI) * std::sin(y * 2 * M_PI);
  };

  for (unsigned int i = 1; i < M + 1; ++i) {
    for (unsigned int j = 1; j < M + 1; ++j) {
      phi(i, j) = f(i * h, j * h);
      if (i <= midpoint && j >= midpoint) {
        phi(i, j) = 0;
      }
    }
  }

  GridFunction mu = Eigen::MatrixXd::Zero(M + 2, M + 2);
  unsigned iter = 0;
  unsigned max_iter = 20;
  unsigned L0 = 2;
  do {
    mu = solver.multigridIteration(mu, h * h * phi, L0);
    iter++;
  } while (solver.residual(l, mu, h * h * phi).cwiseAbs().maxCoeff() > 1e-9 &&
           iter < max_iter);
  ASSERT_NEAR(solver.residual(l, mu, h * h * phi).cwiseAbs().maxCoeff(), 0,
              1e-9);
  ASSERT_NEAR((mu - phi / ((8 * M_PI * M_PI + 1))).cwiseAbs().maxCoeff(), 0,
              1e-4);
}

}  // namespace MFMGLshape::test