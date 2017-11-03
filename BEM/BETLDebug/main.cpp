//#define EIGEN_USE_MKL_ALL
//#include <Eigen/Core>
//#include <Eigen/Dense>

// system includes -------------------------------------------------------------
#include <iostream>
#include <memory>
#include <fstream>
#include <cmath>
#include <complex>

// eth includes ----------------------------------------------------------------
#include <input_interface/input_interface.hpp>
#include <grid_utils/grid_view_factory.hpp>
#include <eth_base/ref_el_types_i.hpp>
#include <eth_base/timer.hpp>

// betl2 includes --------------------------------------------------------------
#include <grid/surface_grid.hpp>
#include <grid/grid_view.hpp>
#include <gmsh_input/gmsh_input.hpp>
#include <cmdl_parser/cmdl_parser.hpp>

// quadrature related includes -------------------------------------------------
#include <fe/febasis.hpp>
#include <fe/fe_enumerators.hpp>
#include <fe/dof_handler.hpp>
#include <bem_integration/integration_traits.hpp>
#include <bem_integration/galerkin_singularity_detector.hpp>
#include <bem_integration/galerkin_quadrature_rules.hpp>
#include <bem_integration/cache.hpp>
#include <bem_integration/galerkin_kernel.hpp>
#include <bem_integration/galerkin_integrator.hpp>

// fundamental solutions -------------------------------------------------------
#include <fundsol/fundsol.hpp>
#include <analytical_functions/fundsol_functor.hpp>

// bem operator ----------------------------------------------------------------
#include <bem_operator/bem_operator.hpp>

// sparse operators ------------------------------------------------------------
#include <sparse_operators/identity_operator.hpp>
#include <sparse_operators/combinatorial_divergence.hpp>
#include <sparse_operators/combinatorial_gradient.hpp>

// utility functions -----------------------------------------------------------
#include <utils/make_matrix.hpp>

// grid function related includes ----------------------------------------------
#include <functional/analytical_grid_function.hpp>
#include <functional/interpolation_grid_function.hpp>
#include <functional/grid_function_operations.hpp>
#include <functional/dof_interpolator.hpp>

// export datasets to vtu-files ------------------------------------------------
#include <vtu_exporter/vtu_exporter.hpp>

// compute h -------------------------------------------------------------------
#include <utils/element_metrics.hpp>

using namespace betl2;
namespace big = betl2::input::gmsh;

template< int ORDER, typename QUADRATURE_RULE_T >
struct Simulate
{
  template< typename GRID_FACTORY_T >
  void operator()( const GRID_FACTORY_T& gridFactory, std::string basename ) const;
};

//------------------------------------------------------------------------------
int main( int argc, char* argv[] )
{
  // parse the command line
  const std::string basename = betl2::parseCommandLine( argc, argv );

  // create input
  std::cout << "Input from: " << basename << ".msh" << std::endl;
  big::Input input( basename );
  std::cout << input << std::endl;
    
  // wrap input interface around the given input
  typedef betl2::input::InputInterface< big::Input > inpInterface_t;
  inpInterface_t inpInterface( input );
    
  // instantiate the grid implementation
  typedef betl2::surfaceGrid::hybrid::Grid grid_t;
  typedef std::shared_ptr< eth::grid::Grid<grid_t::gridTraits_t> > grid_ptr_t;
  grid_ptr_t grid_ptr( new grid_t( inpInterface ) );
  
  // create a grid view factory
  typedef eth::grids::utils::GridViewFactory< grid_t, eth::grid::GridViewTypes::LeafView > grid_factory_t;
  grid_factory_t grid_factory( grid_ptr );

  // commpute global mesh size
  const double h = utils::globalMeshSize( grid_factory );
  std::cout << "Global mesh size h: " << h << std::endl;

  // define the quadrature rule
  const int numQuadPointsTria = 12;
  const int numQuadPointsQuad = 12;
  const int numSingQuadPoints = 12;

  std::cout << "No. of Gaussian points (regular)  = " << numQuadPointsQuad << std::endl
            << "No. of Gaussian points (singular) = " << numSingQuadPoints << std::endl;

  typedef bem::GalerkinQuadratureRule<numQuadPointsTria,numQuadPointsQuad,numSingQuadPoints> quadrature_rule_t;

  {
    const fe::ApproxOrder order = fe::Linear;
    std::cout << "--> APPROX. ORDER = " << order << " <--" << std::endl;
    Simulate< order, quadrature_rule_t >()( grid_factory, basename );
  }

  {
    const fe::ApproxOrder order = fe::Quadratic;
    std::cout << "--> APPROX. ORDER = " << order << " <--" << std::endl;
    Simulate< order, quadrature_rule_t >()( grid_factory, basename );
  }

  //const fe::ApproxOrder order = fe::Cubic;
  //Simulate< order, quadrature_rule_t >()( grid_factory, basename );

  // that's it!
  return EXIT_SUCCESS;
}



template< int ORDER, typename QUADRATURE_RULE_T >
template< typename GRID_FACTORY_T >
void Simulate< ORDER, QUADRATURE_RULE_T >::operator()( const GRID_FACTORY_T& gridFactory, std::string basename ) const
{
  //============================================================================
  //
  // THE DISCRETE FE-SPACES
  //
  //============================================================================
  
  // some abbreviations for convenience
  typedef GRID_FACTORY_T grid_factory_t;

  // define bases
  typedef fe::FEBasis< ORDER-1, fe::FEBasisType::Lagrange > feb_lagrange0_t;
  typedef fe::FEBasis< ORDER  , fe::FEBasisType::Lagrange > feb_lagrange1_t;
  typedef fe::FEBasis< ORDER  , fe::FEBasisType::Div      > feb_div_t;
  
  // define the dofhandler types
  typedef betl2::fe::DofHandler< feb_lagrange0_t,fe::FESContinuity::Discontinuous,
                                 grid_factory_t > dh_lagrange0_t;

  typedef betl2::fe::DofHandler< feb_lagrange1_t,fe::FESContinuity::Continuous,
                                 grid_factory_t > dh_lagrange1_t;

  typedef betl2::fe::DofHandler< feb_div_t,fe::FESContinuity::Continuous,
                                 grid_factory_t > dh_div_t;

  // instantiate dofhandler objects
  dh_lagrange0_t dh_lagrange0;
  dh_lagrange1_t dh_lagrange1;
  dh_div_t       dh_div;

  // distribute the degrees of freedom
  dh_lagrange0.distributeDofs( gridFactory );
  dh_lagrange1.distributeDofs( gridFactory );
  dh_div.distributeDofs      ( gridFactory );

  const int num_dofs_lagrange0 = dh_lagrange0.numDofs( );
  const int num_dofs_lagrange1 = dh_lagrange1.numDofs( );
  const int num_dofs_div       = dh_div.numDofs( );

  std::cout << "Number of created dofs(lagrange1)  = " << num_dofs_lagrange1  << std::endl;
  std::cout << "Number of created dofs(lagrange0)  = " << num_dofs_lagrange0  << std::endl;
  std::cout << "Number of created dofs(div)        = " << num_dofs_div        << std::endl;


  //============================================================================
  //
  // THE RUNTIME CACHE
  //
  //============================================================================
  // define the runtime cache
  constexpr int numQuadPointsTria = QUADRATURE_RULE_T::triaTriaRule_t::getOuterRE( );
  constexpr int numQuadPointsQuad = QUADRATURE_RULE_T::quadQuadRule_t::getOuterRE( );
  typedef bem::Cache< grid_factory_t, numQuadPointsTria, numQuadPointsQuad > cache_t;
  typedef eth::base::Timer< eth::base::TimeUnit::sec > sec_timer_t;
  sec_timer_t timer_cache;
  cache_t cache( gridFactory );
  std::cout << "Cache instantiation took: " << timer_cache << std::endl;


  //============================================================================
  //
  // FUNDAMENTAL SOLUTIONS
  //
  //============================================================================
  typedef bem::FundSol< bem::FSType::LAPLACE, 3 > laplace_fs_t;
  laplace_fs_t laplace_fs;

  //============================================================================
  //
  // GALERKIN KERNELS
  //
  //============================================================================
  typedef bem::GalerkinKernel< laplace_fs_t, bem::FSLayer::SL, feb_lagrange0_t, feb_lagrange0_t > lagr_sl_kernel_t;
  typedef bem::GalerkinKernel< laplace_fs_t, bem::FSLayer::DL, feb_lagrange0_t, feb_lagrange1_t > lagr_dl_kernel_t;
  typedef bem::GalerkinKernel< laplace_fs_t, bem::FSLayer::SL, feb_div_t      , feb_div_t       > lagr_div_kernel_t; 
  lagr_sl_kernel_t  lagr_sl_kernel ( laplace_fs );
  lagr_dl_kernel_t  lagr_dl_kernel ( laplace_fs );
  lagr_div_kernel_t lagr_div_kernel( laplace_fs );
  
  //============================================================================
  //
  // THE SINGULARITY DETECTOR
  //
  //============================================================================
  typedef bem::GalerkinSingularityDetector<grid_factory_t> singularity_detector_t;
  singularity_detector_t singularity_detector( gridFactory );


  //============================================================================
  //
  // BEM INTEGRATORS
  //
  //============================================================================

  // The following boolean value is optional. If not given it is assumed to be true.
  // It indicates whether the gram determinants are to be included into the integral or not.
  // In general, these determinants are mandatory. However, certain linear combinations
  // of scalar potentials might be used to express potentials related to edge functions.
  // In this particular case, the gram determinants cancel out and should(must!) be neglected.
  const bool withGram = true; 
  typedef bem::IntegrationTraits< QUADRATURE_RULE_T,
                                  grid_factory_t, 
                                  singularity_detector_t, 
                                  cache_t, 
                                  withGram, /* on element X */
                                  withGram  /* on element Y */ > integrationTraits;


  typedef bem::GalerkinIntegrator< lagr_sl_kernel_t,integrationTraits > lagr_sl_integrator_t;
  typedef bem::GalerkinIntegrator< lagr_dl_kernel_t,integrationTraits > lagr_dl_integrator_t;
  // an integrator for edge functions. In this case the withGram-flag within 'integrationTraits' 
  // has no effect since integrators for edge-functions neglect the Gram determinants in any case
  typedef bem::GalerkinIntegrator< lagr_div_kernel_t,integrationTraits > lagr_div_integrator_t;

  sec_timer_t timer_int;
  lagr_sl_integrator_t  lagr_sl_integrator ( lagr_sl_kernel , singularity_detector, cache );
  lagr_dl_integrator_t  lagr_dl_integrator ( lagr_dl_kernel , singularity_detector, cache );
  lagr_div_integrator_t lagr_div_integrator( lagr_div_kernel, singularity_detector, cache );
  std::cout << "Integrator instantiation took: " << timer_int << std::endl;



  //============================================================================
  //
  // BEM OPERATORS
  //
  //============================================================================
  sec_timer_t timer_bem;
  typedef BemOperator< lagr_sl_integrator_t,
                       typename dh_lagrange0_t::fespace_t > bem_op_V_t;
  typedef BemOperator< lagr_dl_integrator_t,
                       typename dh_lagrange0_t::fespace_t,
                       typename dh_lagrange1_t::fespace_t > bem_op_K_t;
  typedef BemOperator< lagr_div_integrator_t,
                       typename dh_div_t::fespace_t > bem_op_div_t;

  bem_op_V_t   bem_op_V  ( lagr_sl_integrator , dh_lagrange0.fespace() );
  bem_op_K_t   bem_op_K  ( lagr_dl_integrator , dh_lagrange0.fespace(), dh_lagrange1.fespace() );
  bem_op_div_t bem_op_div( lagr_div_integrator, dh_div.fespace() );
  
  sec_timer_t timer_op_V;
  bem_op_V.compute( );
  std::cout << "Computation of <V phi, phi>, phi in H^{-1/2} took: " << timer_op_V << std::endl;
  sec_timer_t timer_op_K;
  bem_op_K.compute( );
  std::cout << "Computation of <K psi, phi>, phi in H^{-1/2}, psi in H^{1/2} took: " 
            << timer_op_V << std::endl;
  sec_timer_t timer_op_div;
  bem_op_div.compute( );
  std::cout << "Computation of <V u, u>, u in H^{-1/2}(div) took: " << timer_op_div << std::endl;
  std::cout << "Computation of bem operators took: " << timer_bem << std::endl;


  //============================================================================
  //
  // CREATE MASS-MATRICES
  //
  //============================================================================
  typedef IdentityOperator< typename dh_lagrange0_t::fespace_t,
                            typename dh_lagrange1_t::fespace_t > discrete_op_t;
  discrete_op_t discrete_op( dh_lagrange0.fespace(), dh_lagrange1.fespace() );
  discrete_op.compute( );
  const auto& M = discrete_op.matrix();

  
  //============================================================================
  //
  // THE DISCRETE EMBEDDING H(lagrange) -> H(div)
  //
  //============================================================================
  typedef CombinatorialGradient< typename dh_lagrange1_t::fespace_t,
                                 typename dh_div_t      ::fespace_t > combinatorial_gradient_t;
  combinatorial_gradient_t combinatorial_gradient( dh_lagrange1.fespace(), dh_div.fespace() );
  combinatorial_gradient.compute( cache );

  const auto& C = combinatorial_gradient.matrix();



  //============================================================================
  //
  // CREATE ANALYTICAL GRID FUNCTION
  //
  //============================================================================
  typedef analytical::FundsolFunctor< laplace_fs_t > fundsol_functor_t;

  typedef bem::AnalyticalGridFunction< grid_factory_t, fundsol_functor_t, Trace::Dirichlet > analytical_dirichlet_t;
  typedef bem::AnalyticalGridFunction< grid_factory_t, fundsol_functor_t, Trace::Neumann   > analytical_neumann_t;

  typedef utils::MakeMatrix<double,3,1> matrix_maker;
  const auto source = matrix_maker()( { 1.1, 1.2, 1.03 } );
  
  const fundsol_functor_t      fundsol_functor( laplace_fs, source );
  const analytical_dirichlet_t analytical_dirichlet( gridFactory, fundsol_functor );
  const analytical_neumann_t   analytical_neumann  ( gridFactory, fundsol_functor );


  //============================================================================
  //
  // CREATE COEFFICIENTS FROM GRID-FUNCTIONS
  //
  //============================================================================
  const auto& fes_lagrange1 = dh_lagrange1.fespace( );
  const auto  coeff_gD      = DofInterpolator()( analytical_dirichlet, fes_lagrange1 );

  const auto& fes_lagrange0 = dh_lagrange0.fespace( );
  const auto  coeff_gN      = DofInterpolator()( analytical_neumann, fes_lagrange0 );

  //============================================================================
  //
  // CREATE INTERPOLATING FUNCTIONS BASED ON THE DISCRETE DOFS
  //
  //============================================================================
  typedef InterpolationGridFunction< grid_factory_t,typename dh_lagrange1_t::fespace_t,double> interpolation_dirichlet_t;
  typedef InterpolationGridFunction< grid_factory_t,typename dh_lagrange0_t::fespace_t,double> interpolation_neumann_t;


  interpolation_dirichlet_t interpolation_dirichlet( gridFactory, dh_lagrange1.fespace(), coeff_gD );
  interpolation_neumann_t   interpolation_neumann  ( gridFactory, dh_lagrange0.fespace(), coeff_gN );


  //============================================================================
  //
  // COMPUTE RESIDUALS
  //
  //============================================================================
  
  // lhs := V u_N - K u_D
  const auto& V    = bem_op_V.matrix();
  const auto& K    = bem_op_K.matrix();
  const auto& Vdiv = bem_op_div.matrix();

  typedef Eigen::Matrix<double,Eigen::Dynamic,1> vector_t;
  const vector_t lhs_D = V * coeff_gN - K * coeff_gD;
  const vector_t lhs_N = K.transpose() * coeff_gN + C * Vdiv * C.transpose() * coeff_gD;

  // rhs := 0.5 M u_D
  const vector_t rhs_D = 0.5 * M * coeff_gD;
  const vector_t rhs_N = 0.5 * M.transpose() * coeff_gN;

  // output Euclidean norm |lhs - rhs|
  std::cout << "Euclidean norm of residuum: | V uN - K uD - 0.5 M uD | = " << ( lhs_D - rhs_D ).norm( ) << std::endl;
  std::cout << "Euclidean norm of residuum: |K' uN - 0.5 M uN + W uD | = " << ( lhs_N - rhs_N ).norm( ) << std::endl;


  //============================================================================
  //
  // CHECK ANALYTICAL VS. INTERPOLATED GRID-FUNCTIONS
  //
  //============================================================================
  const double UD_nrm = analytical_dirichlet.norm( );
  const double UN_nrm = analytical_neumann.norm( );
  const double UD_diff_nrm = ( analytical_dirichlet - interpolation_dirichlet ).norm( );
  const double UN_diff_nrm = ( analytical_neumann   - interpolation_neumann   ).norm( );
  
  std::cout << "           || UD - UD_h ||_2 = " << UD_diff_nrm << std::endl
            << "           || UN - UN_h ||_2 = " << UN_diff_nrm << std::endl
            << "|| UD - UD_h ||_2 / ||UD||_2 = " << UD_diff_nrm / UD_nrm << std::endl
            << "|| UN - UN_h ||_2 / ||UN||_2 = " << UN_diff_nrm / UN_nrm << std::endl;


  //============================================================================
  //
  // EXPORTER SECTION
  //
  //============================================================================
  typedef vtu::Exporter< grid_factory_t > exporter_t;
  exporter_t exporter( gridFactory, basename );
  
  exporter
    ( "UD   [exact]"    , analytical_dirichlet   , vtu::Entity::Point )
    ( "UN   [exact]"    , analytical_neumann     , vtu::Entity::Cell  )
    ( "UD_h [interpl'd]", interpolation_dirichlet, vtu::Entity::Point )
    ( "UN_h [interpl'd]", interpolation_neumann  , vtu::Entity::Cell  );
}
