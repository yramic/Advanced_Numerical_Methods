////
//// Copyright (C) 2017 SAM (D-MATH) @ ETH Zurich
//// Author(s): curzuato < >
//// Contributors:  dcasati
//// This file is part of the AdvNumCSE repository.
////
#include <Eigen/Core>
#include <Eigen/Dense>
// system includes -------------------------------------------------------------
#include <cmath>
#include <fstream>
#include <iostream>
#include <memory>
// eth includes ----------------------------------------------------------------
#include <eth_base/ref_el_types_i.hpp>
#include <eth_base/timer.hpp>
#include <grid_utils/grid_view_factory.hpp>
#include <input_interface/input_interface.hpp>
// betl2 includes --------------------------------------------------------------
#include <cmdl_parser/cmdl_parser.hpp>
#include <gmsh_input/gmsh_input.hpp>
#include <grid/grid_view.hpp>
#include <grid/surface_grid.hpp>
// quadrature related includes -------------------------------------------------
#include <bem_integration/cache.hpp>
#include <bem_integration/galerkin_integrator.hpp>
#include <bem_integration/galerkin_kernel.hpp>
#include <bem_integration/galerkin_quadrature_rules.hpp>
#include <bem_integration/galerkin_singularity_detector.hpp>
#include <bem_integration/integration_traits.hpp>
#include <fe/dof_handler.hpp>
#include <fe/fe_enumerators.hpp>
#include <fe/febasis.hpp>
// fundamental solutions -------------------------------------------------------
#include <analytical_functions/fundsol_functor.hpp>
#include <fundsol/fundsol.hpp>
// bem operator ----------------------------------------------------------------
#include <bem_operator/bem_operator.hpp>
// sparse operators ------------------------------------------------------------
#include <sparse_operators/combinatorial_gradient.hpp>
#include <sparse_operators/identity_operator.hpp>
// utility functions -----------------------------------------------------------
#include <utils/make_matrix.hpp>
// grid function related includes ----------------------------------------------
#include <functional/analytical_grid_function.hpp>
#include <functional/dof_interpolator.hpp>
#include <functional/grid_function_operations.hpp>
#include <functional/interpolation_grid_function.hpp>

using namespace betl2;
namespace big = betl2::input::gmsh;

// Define types to be used
// Grid
typedef betl2::input::InputInterface<big::Input> inpInterface_t;
typedef betl2::surfaceGrid::hybrid::Grid grid_t;
typedef std::shared_ptr<eth::grid::Grid<grid_t::gridTraits_t> > grid_ptr_t;
typedef eth::grids::utils::GridViewFactory<grid_t,
                                           eth::grid::GridViewTypes::LeafView>
    grid_factory_t;
// FE-Basis
const fe::ApproxOrder order = fe::Linear;
typedef fe::FEBasis<order - 1, fe::FEBasisType::Lagrange> feb_lagr0_t;
typedef fe::FEBasis<order, fe::FEBasisType::Lagrange> feb_lagr1_t;
typedef fe::FEBasis<order, fe::FEBasisType::Div> feb_div_t;
// Dof-handlers
typedef betl2::fe::DofHandler<feb_lagr0_t, fe::FESContinuity::Discontinuous,
                              grid_factory_t>
    dh_lagrange0_t;
typedef betl2::fe::DofHandler<feb_lagr1_t, fe::FESContinuity::Continuous,
                              grid_factory_t>
    dh_lagrange1_t;
typedef betl2::fe::DofHandler<feb_div_t, fe::FESContinuity::Continuous,
                              grid_factory_t>
    dh_div_t;
// Fundamental solution
typedef bem::FundSol<bem::FSType::LAPLACE, 3> laplace_fs_t;
// Singularity detector
typedef bem::GalerkinSingularityDetector<grid_factory_t> singularity_detector_t;
// Eigen
typedef Eigen::MatrixXd matrix_t;

//------------------------------------------------------------------------------
void computeV(const grid_factory_t& gridFactory,
              const dh_lagrange0_t& dh_lagrange0,
              const laplace_fs_t& laplace_fs, matrix_t& V) {
  // TODO: Compute V
}

//------------------------------------------------------------------------------
void computeW(const grid_factory_t& gridFactory,
              const dh_lagrange1_t& dh_lagrange1, const dh_div_t& dh_div,
              const laplace_fs_t& laplace_fs, matrix_t& W) {
  // TODO: Compute W
}

//------------------------------------------------------------------------------
void computeK(const grid_factory_t& gridFactory,
              const dh_lagrange0_t& dh_lagrange0,
              const dh_lagrange1_t& dh_lagrange1,
              const laplace_fs_t& laplace_fs, matrix_t& K) {
  // GALERKIN KERNEL
  typedef bem::GalerkinKernel<laplace_fs_t, bem::FSLayer::DL, feb_lagr0_t,
                              feb_lagr1_t>
      lagr_dl_kernel_t;
  lagr_dl_kernel_t lagr_dl_kernel(laplace_fs);

  // THE SINGULARITY DETECTOR
  singularity_detector_t singularity_detector(gridFactory);

  // DEFINE QUADRATURE
  const int numQPtsTria = 12;
  const int numQPtsQuad = 12;
  const int numSingQPts = 12;

  typedef bem::GalerkinQuadratureRule<numQPtsTria, numQPtsQuad, numSingQPts>
      quadrature_rule_t;

  // THE RUNTIME CACHE
  typedef bem::Cache<grid_factory_t, numQPtsTria, numQPtsQuad> cache_t;
  cache_t cache(gridFactory);

  // BEM INTEGRATOR
  const bool withGram = true;
  typedef bem::IntegrationTraits<quadrature_rule_t, grid_factory_t,
                                 singularity_detector_t, cache_t,
                                 withGram, /* on element X */
                                 withGram /* on element Y */>
      integrationTraits;

  typedef bem::GalerkinIntegrator<lagr_dl_kernel_t, integrationTraits>
      lagr_dl_integrator_t;
  lagr_dl_integrator_t lagr_dl_integrator(lagr_dl_kernel, singularity_detector,
                                          cache);

  // BEM OPERATOR
  typedef BemOperator<lagr_dl_integrator_t, typename dh_lagrange0_t::fespace_t,
                      typename dh_lagrange1_t::fespace_t>
      bem_op_K_t;
  bem_op_K_t bem_op_K(lagr_dl_integrator, dh_lagrange0.fespace(),
                      dh_lagrange1.fespace());
  bem_op_K.compute();
  K = bem_op_K.matrix();
}

//------------------------------------------------------------------------------
void debugV(const matrix_t& V) {
  // Check V is spd by means of its eigenvalues

  // TODO: Implement your code
}

//------------------------------------------------------------------------------
void debugW(const matrix_t& W) {
  // TODO: Implement your code
}

//------------------------------------------------------------------------------
Eigen::VectorXd computeDirichletResidual(const grid_factory_t& gridFactory,
                                         const dh_lagrange0_t& dh_lagrange0,
                                         const dh_lagrange1_t& dh_lagrange1,
                                         const laplace_fs_t& laplace_fs) {
  // TODO: Implement your code

  // CREATE ANALYTICAL GRID FUNCTION FOR TRACES
  typedef analytical::FundsolFunctor<laplace_fs_t> fundsol_functor_t;

  typedef bem::AnalyticalGridFunction<grid_factory_t, fundsol_functor_t,
                                      Trace::Dirichlet>
      analytical_dirichlet_t;
  typedef bem::AnalyticalGridFunction<grid_factory_t, fundsol_functor_t,
                                      Trace::Neumann>
      analytical_neumann_t;

  typedef utils::MakeMatrix<double, 3, 1> matrix_maker;
  const auto source = matrix_maker()({1.1, 1.2, 1.03});

  const fundsol_functor_t fundsol_functor(laplace_fs, source);
  const analytical_dirichlet_t analytical_dirichlet(gridFactory,
                                                    fundsol_functor);
  const analytical_neumann_t analytical_neumann(gridFactory, fundsol_functor);

  // CREATE COEFFICIENTS VECTOR FOR GRID-FUNCTIONS OF THE TRACES
  const auto coeff_gD =
      DofInterpolator()(analytical_dirichlet, dh_lagrange1.fespace());

  const auto coeff_gN =
      DofInterpolator()(analytical_neumann, dh_lagrange0.fespace());

  // TODO: Implement your code
  Eigen::VectorXd res_D;

  return res_D;
}

//------------------------------------------------------------------------------
Eigen::VectorXd computeNeumannResidual(const grid_factory_t& gridFactory,
                                       const dh_lagrange0_t& dh_lagrange0,
                                       const dh_lagrange1_t& dh_lagrange1,
                                       const dh_div_t& dh_div,
                                       const laplace_fs_t& laplace_fs) {
  // TODO: Implement your code

  // CREATE ANALYTICAL GRID FUNCTION FOR TRACES
  typedef analytical::FundsolFunctor<laplace_fs_t> fundsol_functor_t;

  typedef bem::AnalyticalGridFunction<grid_factory_t, fundsol_functor_t,
                                      Trace::Dirichlet>
      analytical_dirichlet_t;
  typedef bem::AnalyticalGridFunction<grid_factory_t, fundsol_functor_t,
                                      Trace::Neumann>
      analytical_neumann_t;

  typedef utils::MakeMatrix<double, 3, 1> matrix_maker;
  const auto source = matrix_maker()({1.1, 1.2, 1.03});

  const fundsol_functor_t fundsol_functor(laplace_fs, source);
  const analytical_dirichlet_t analytical_dirichlet(gridFactory,
                                                    fundsol_functor);
  const analytical_neumann_t analytical_neumann(gridFactory, fundsol_functor);

  // CREATE COEFFICIENTS VECTOR FOR GRID-FUNCTIONS OF THE TRACES
  const auto coeff_gD =
      DofInterpolator()(analytical_dirichlet, dh_lagrange1.fespace());

  const auto coeff_gN =
      DofInterpolator()(analytical_neumann, dh_lagrange0.fespace());

  // TODO: Implement your code
  Eigen::VectorXd res_N;

  return res_N;
}

//------------------------------------------------------------------------------
int main(int argc, char* argv[]) {
  Eigen::VectorXi levels(4);
  Eigen::VectorXd rDNorm(4), rNNorm(4);
  levels << 32, 128, 512, 2048;
  const std::string path = "../BEM/BETL-Debug/meshes/";

  for (int k = 0; k < 4; k++) {
    //============================================================================
    // READ MESH
    //============================================================================
    // Create input from mesh

    const std::string basename = "sphere_" + std::to_string(levels(k));
    std::cout << "Input from: " << basename << ".msh" << std::endl;
    big::Input input(path + basename);

    // Wrap input interface around the given input
    inpInterface_t inpInterface(input);

    //============================================================================
    // CREATE GRID FROM INPUT MESH
    //============================================================================
    // Instantiate the grid implementation
    grid_ptr_t grid_ptr(new grid_t(inpInterface));

    // Create a grid view factory
    grid_factory_t gridFactory(grid_ptr);

    //============================================================================
    // THE DISCRETE FE-SPACES
    //============================================================================
    // - Instantiate dofhandler objects
    dh_lagrange0_t dh_lagrange0;
    dh_lagrange1_t dh_lagrange1;
    dh_div_t dh_div;

    // - Distribute the degrees of freedom
    dh_lagrange0.distributeDofs(gridFactory);
    dh_lagrange1.distributeDofs(gridFactory);
    dh_div.distributeDofs(gridFactory);
    std::cout << "Created " << dh_lagrange1.numDofs() << " dofs (lagrange1).\n"
              << "Created " << dh_lagrange0.numDofs() << " dofs (lagrange0).\n"
              << "Created " << dh_div.numDofs() << " dofs (div)." << std::endl;

    //============================================================================
    // FUNDAMENTAL SOLUTION FOR LAPLACE EQUATION
    //============================================================================
    laplace_fs_t laplace_fs;

    //============================================================================
    // TEST FOR V AND W
    //============================================================================
    matrix_t V, W;
    computeV(gridFactory, dh_lagrange0, laplace_fs, V);
    debugV(V);

    computeW(gridFactory, dh_lagrange1, dh_div, laplace_fs, W);
    debugW(W);

    //============================================================================
    // RESIDUALS
    //============================================================================

    // TODO: Implement your code
  }

  // Output
  std::ofstream out_rdNorm("BETL-Debug_rDnorm.txt");
  out_rdNorm << rDNorm;
  out_rdNorm.close();

  std::ofstream out_rNNorm("BETL-Debug_rNnorm.txt");
  out_rNNorm << rNNorm;
  out_rNNorm.close();

  std::ofstream out_N("BETL-Debug_levels.txt");
  out_N << levels;
  out_N.close();

  // that's it!
  return EXIT_SUCCESS;
}
