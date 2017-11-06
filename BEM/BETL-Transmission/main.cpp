#include <Eigen/Core>
#include <Eigen/Dense>
// system includes -------------------------------------------------------------
#include <iostream>
#include <memory>
#include <fstream>
#include <cmath>
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
#include <sparse_operators/combinatorial_gradient.hpp>
// utility functions -----------------------------------------------------------
#include <utils/make_matrix.hpp>
// grid function related includes ----------------------------------------------
#include <functional/analytical_grid_function.hpp>
#include <functional/interpolation_grid_function.hpp>
#include <functional/grid_function_operations.hpp>
#include <functional/dof_interpolator.hpp>
//#include <functional/L_two_product_evaluator.hpp>
// own includes ----------------------------------------------------------------
#include "transmission_system_matrix.hpp"


using namespace betl2;
namespace big = betl2::input::gmsh;

//------------------------------------------------------------------------------
void solveTransmissionProblemSphere(const double& alpha, const int& N){
  //============================================================================
  // READ MESH
  //============================================================================
  // Create input from mesh
  const std::string path = "../BEM/BETL-Transmission/meshes/";
  const std::string basename ="sphere_" + std::to_string(N);
  std::cout << "Input from: " << basename << ".msh" << std::endl;
  big::Input input( path+basename );
  
  // Wrap input interface around the given input
  typedef betl2::input::InputInterface< big::Input > inpInterface_t;
  inpInterface_t inpInterface( input );

  
  //============================================================================
  // CREATE GRID FROM INPUT MESH
  //============================================================================
  // Instantiate the grid implementation
  typedef betl2::surfaceGrid::hybrid::Grid grid_t;
  typedef std::shared_ptr< eth::grid::Grid<grid_t::gridTraits_t> > grid_ptr_t;
  grid_ptr_t grid_ptr( new grid_t( inpInterface ) );
  
  // Create a grid view factory
  typedef eth::grids::utils::GridViewFactory< grid_t,
					      eth::grid::GridViewTypes::LeafView
					      > grid_factory_t;
  grid_factory_t gridFactory( grid_ptr );

  
  //============================================================================
  // THE DISCRETE FE-SPACES
  //============================================================================
  const fe::ApproxOrder order = fe::Linear;
  typedef fe::FEBasis< order-1, fe::FEBasisType::Lagrange > feb_lagr0_t;
  typedef fe::FEBasis< order  , fe::FEBasisType::Lagrange > feb_lagr1_t;
  typedef fe::FEBasis< order  , fe::FEBasisType::Div      > feb_div_t;
  
  // 2 - Define the dofhandler types
  typedef betl2::fe::DofHandler< feb_lagr0_t,fe::FESContinuity::Discontinuous,
                                 grid_factory_t > dh_lagrange0_t;
  typedef betl2::fe::DofHandler< feb_lagr1_t,fe::FESContinuity::Continuous,
                                 grid_factory_t > dh_lagrange1_t;
  typedef betl2::fe::DofHandler< feb_div_t,fe::FESContinuity::Continuous,
				 grid_factory_t > dh_div_t;
  
  // - Instantiate dofhandler objects
  dh_lagrange0_t dh_lagrange0;
  dh_lagrange1_t dh_lagrange1;
  dh_div_t       dh_div;

  // - Distribute the degrees of freedom
  dh_lagrange0.distributeDofs( gridFactory );
  dh_lagrange1.distributeDofs( gridFactory );
  dh_div.distributeDofs      ( gridFactory );
  std::cout << "Created " << dh_lagrange1.numDofs() << " dofs (lagrange1).\n"
	    << "Created " << dh_lagrange0.numDofs() << " dofs (lagrange0).\n"
	    << "Created " << dh_div.numDofs() << " dofs (div)." << std::endl;


  //============================================================================
  // FUNDAMENTAL SOLUTION FOR LAPLACE EQUATION
  //============================================================================
  typedef bem::FundSol< bem::FSType::LAPLACE, 3 > laplace_fs_t;
  laplace_fs_t laplace_fs;

  
  //============================================================================
  // GALERKIN KERNELS
  //============================================================================
  // - Kernel for V: Single layer with $\qbe$ for test and trial spaces
  typedef bem::GalerkinKernel< laplace_fs_t, bem::FSLayer::SL,
			       feb_lagr0_t, feb_lagr0_t > lagr_sl_kernel_t;
  // - Kernel for K: Double layer with $\qbe, \sbe$ bases (we get K' by transposing)
  typedef bem::GalerkinKernel< laplace_fs_t, bem::FSLayer::DL,
			       feb_lagr0_t, feb_lagr1_t > lagr_dl_kernel_t;
  // - Kernel for W using integration by parts
  typedef bem::GalerkinKernel< laplace_fs_t, bem::FSLayer::SL,
			       feb_div_t  , feb_div_t   > div_sl_kernel_t;  
  lagr_sl_kernel_t lagr_sl_kernel( laplace_fs );
  lagr_dl_kernel_t lagr_dl_kernel( laplace_fs );
  div_sl_kernel_t  div_sl_kernel ( laplace_fs );

  
  //============================================================================
  // THE SINGULARITY DETECTOR
  //============================================================================
  typedef bem::GalerkinSingularityDetector<grid_factory_t> singularity_detector_t;
  singularity_detector_t singularity_detector( gridFactory );


  //============================================================================
  // DEFINE QUADRATURE
  //============================================================================
  const int numQPtsTria = 12;
  const int numQPtsQuad = 12;
  const int numSingQPts = 12;

  typedef bem::GalerkinQuadratureRule<numQPtsTria, numQPtsQuad,
				      numSingQPts> quadrature_rule_t;

  
  //============================================================================
  // THE RUNTIME CACHE
  //============================================================================
  // Define the runtime cache
  typedef bem::Cache< grid_factory_t, numQPtsTria, numQPtsQuad > cache_t;
  typedef eth::base::Timer< eth::base::TimeUnit::sec > sec_timer_t;
  sec_timer_t timer_cache;
  cache_t cache( gridFactory );
  std::cout << "Cache instantiation took: " << timer_cache << std::endl;

  
  //============================================================================
  // BEM INTEGRATORS
  //============================================================================
  // Define integrator to be used with the kernels above.
  typedef bem::IntegrationTraits< quadrature_rule_t,
                                  grid_factory_t, 
                                  singularity_detector_t, 
                                  cache_t > integrationTraits;

  typedef bem::GalerkinIntegrator< lagr_sl_kernel_t,
				   integrationTraits > lagr_sl_integrator_t;
  typedef bem::GalerkinIntegrator< lagr_dl_kernel_t,
				   integrationTraits > lagr_dl_integrator_t;
  typedef bem::GalerkinIntegrator< div_sl_kernel_t,
				   integrationTraits > div_sl_integrator_t;

  sec_timer_t timer_int;
  lagr_sl_integrator_t lagr_sl_integrator( lagr_sl_kernel, singularity_detector, cache );
  lagr_dl_integrator_t lagr_dl_integrator( lagr_dl_kernel, singularity_detector, cache );
  div_sl_integrator_t  div_sl_integrator ( div_sl_kernel , singularity_detector, cache );
  std::cout << "Integrator instantiation took: " << timer_int << std::endl;


  //============================================================================
  // BEM OPERATORS
  //============================================================================
  // - V 
  typedef BemOperator< lagr_sl_integrator_t,
                       typename dh_lagrange0_t::fespace_t > bem_op_V_t;
  bem_op_V_t   bem_op_V  ( lagr_sl_integrator , dh_lagrange0.fespace() );
  bem_op_V.compute( );
  const auto& V    = bem_op_V.matrix();

  // - K
  typedef BemOperator< lagr_dl_integrator_t,
                       typename dh_lagrange0_t::fespace_t,
                       typename dh_lagrange1_t::fespace_t > bem_op_K_t;
  bem_op_K_t   bem_op_K  ( lagr_dl_integrator , dh_lagrange0.fespace(),
			   dh_lagrange1.fespace() );
  bem_op_K.compute( );
  const auto& K    = bem_op_K.matrix();

  // - W (using integration by parts):
  typedef BemOperator< div_sl_integrator_t,
                       typename dh_div_t::fespace_t > bem_op_div_t;
  bem_op_div_t bem_op_div( div_sl_integrator, dh_div.fespace() );
  bem_op_div.compute( );
  const auto& Vdiv = bem_op_div.matrix();

  typedef CombinatorialGradient< typename dh_lagrange1_t::fespace_t,
                                 typename dh_div_t      ::fespace_t > cgrad_op_t;
  cgrad_op_t cgrad_op( dh_lagrange1.fespace(), dh_div.fespace() );
  cgrad_op.compute( cache );
  const auto& C = cgrad_op.matrix();

  const auto& W = C * Vdiv * C.transpose();


  //============================================================================
  // SYSTEM MATRIX
  //============================================================================
  typedef TransmissionSystemMatrix TransmissionSystemMatrix_t;
  TransmissionSystemMatrix_t systemMatrix( V, K, W, alpha );
  systemMatrix.compute();
  auto A = systemMatrix.matrix();


  //============================================================================
  // CREATE ANALYTICAL GRID FUNCTION FOR TRACES OF INCIDENT STATIC WAVE
  //============================================================================
  typedef analytical::FundsolFunctor< laplace_fs_t > fundsol_functor_t;

  typedef bem::AnalyticalGridFunction< grid_factory_t, fundsol_functor_t,
				       Trace::Dirichlet > analytical_dirichlet_t;
  typedef bem::AnalyticalGridFunction< grid_factory_t, fundsol_functor_t,
				       Trace::Neumann   > analytical_neumann_t;

  typedef utils::MakeMatrix<double,3,1> matrix_maker;
  const auto source = matrix_maker()( { 0.0, 0.0, 1.0 } );
  
  const fundsol_functor_t      fundsol_functor( laplace_fs, source );
  const analytical_dirichlet_t analytical_dirichlet( gridFactory, fundsol_functor );
  const analytical_neumann_t   analytical_neumann  ( gridFactory, fundsol_functor );


  //============================================================================
  // CREATE COEFFICIENTS VECTOR FOR GRID-FUNCTIONS OF THE TRACES
  //============================================================================
  const auto  coeff_TDui  = DofInterpolator()( analytical_dirichlet,
					       dh_lagrange1.fespace() );

  const auto  coeff_TNui  = DofInterpolator()( analytical_neumann,
					       dh_lagrange0.fespace() );


  //============================================================================
  // BUILD RIGHT HAND SIDE USING TRACES DEFINED ABOVE
  //============================================================================
  // Evaluate <TD uinc, v> and <TN uinc, v>
  /*
  typedef L2ProductEval<grid_factory_t, analytical_dirichlet_t,
			dh_lagrange0_t::fespace_t > L2ProductEvalD_t;
  typedef L2ProductEval<grid_factory_t, analytical_neumann_t,
			dh_lagrange1_t::fespace_t> L2ProductEvalN_t;
  
  const L2ProductEvalD_t& L2ProdEvalDirichlet = L2ProductEvalD_t( gridFactory,
								  analytical_dirichlet,
								  dh_lagrange0.fespace() );
  const L2ProductEvalN_t& L2ProdEvalNeumann = L2ProductEvalN_t( gridFactory,
								analytical_neumann,
							        dh_lagrange1.fespace() );
  const Eigen::VectorXd rhs_D = L2ProdEvalDirichlet.evaluate();
  const Eigen::VectorXd rhs_N = L2ProdEvalNeumann.evaluate();
  */
  // Mass matrices
  typedef IdentityOperator< dh_lagrange0_t::fespace_t,
			    dh_lagrange1_t::fespace_t > id_op01_t;
  id_op01_t id_op01( dh_lagrange0.fespace(), dh_lagrange1.fespace() );
  id_op01.compute( );
  const auto& M = id_op01.matrix();
  // Evaluate <TD uinc, v> and <TN uinc, v>
  const Eigen::VectorXd rhs_D = M* coeff_TDui;
  const Eigen::VectorXd rhs_N = M.transpose()*coeff_TNui;
  Eigen::VectorXd rhs(rhs_D.rows() + rhs_N.rows(),1);
  rhs.setZero();
  rhs.block(0,0,rhs_D.rows(),1) = -rhs_D;
  rhs.block(rhs_D.rows(), 0, rhs_N.rows(), 1) = rhs_N;


  //============================================================================
  // SOLVE (USING DIRECT SOLVER)
  //============================================================================
  // Define the solver
  typedef Eigen::ColPivHouseholderQR< TransmissionSystemMatrix_t::matrix_t
				      > directSolver_t;
  directSolver_t solverA( A );   
  // Initialize the solver
  solverA.compute( A );   
  // Solve the system
  const auto& sol = solverA.solve( rhs );   
  std::cout << "Solved linear system Ax = rhs. " << std::endl;


  //============================================================================
  // OUTPUT
  //============================================================================
  std::ofstream sol_out( "sol_" + basename + ".dat" ); 
  sol_out << std::setprecision(18) << sol;  
  sol_out.close();

}


//------------------------------------------------------------------------------
int main( int argc, char* argv[] )
{
    Eigen::VectorXi levels(4);
    levels << 32, 128, 512, 2048;
    int Nmax = 3;
    double alpha = 2;
    
    for(int k=0; k<Nmax; k++){
      solveTransmissionProblemSphere(alpha, levels(k));
    }

  



  // that's it!
  return EXIT_SUCCESS;
}


