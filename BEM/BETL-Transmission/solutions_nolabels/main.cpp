//// 
//// Copyright (C) 2017 SAM (D-MATH) @ ETH Zurich
//// Author(s): curzuato < > 
//// Contributors:  dcasati 
//// This file is part of the AdvNumCSE repository.
////
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
#include <analytical_functions/linear_excitation_field.hpp>


using namespace betl2;
namespace big = betl2::input::gmsh;

typedef Eigen::MatrixXd matrix_t;
const std::string path = "../../meshes/";


//------------------------------------------------------------------------------
matrix_t TransmissionSystemMatrix( const matrix_t&  V, const matrix_t&  K,
				   const matrix_t&  W, const double & alpha){
  matrix_t A(K.rows() + W.rows(), K.cols() + V.cols());
  A.block(0, 0              , K.rows(), K.cols()) =  2.* K;
  A.block(0, K.cols()       , V.rows(), V.cols()) = -(1./alpha + 1)*V;
  A.block(K.rows(), 0       , W.rows(), W.cols()) = -(alpha + 1.)*W;
  A.block(K.rows(), W.cols(), K.cols(), K.rows()) = - 2.*K.transpose();
  return A;
}


//------------------------------------------------------------------------------
template<typename GRID_FACTORY, typename DH_LAGR0, typename DH_LAGR1>
double computeEnergy(const Eigen::VectorXd& sol, const GRID_FACTORY gridFactory,
		     const DH_LAGR0& dh_lagr0, const DH_LAGR1& dh_lagr1){
  // Extract Dirichlet and Neumann coefficients out of solution vector
  const auto& TDu_h  = sol.segment(0, dh_lagr1.numDofs());
  const auto& TNu_h  = sol.segment(dh_lagr1.numDofs(), dh_lagr0.numDofs());
  
  typedef IdentityOperator<  typename DH_LAGR0::fespace_t,
			     typename DH_LAGR1::fespace_t> id_op01_t;
  id_op01_t id_op01( dh_lagr0.fespace(), dh_lagr1.fespace() );
  id_op01.compute( );
  const auto& M = id_op01.matrix();
  return TNu_h.transpose()*M*TDu_h;
  
}


//------------------------------------------------------------------------------
Eigen::VectorXd solveTransmissionProblem(const std::string meshname,
					 const double& alpha){
  //============================================================================
  // READ MESH
  //============================================================================
  // Create input from mesh
  std::cout << "Input from: " << meshname << ".msh" << std::endl;
  big::Input input( path+meshname );
  
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
  int lagr0_numDofs = dh_lagrange0.numDofs();
  int lagr1_numDofs = dh_lagrange1.numDofs();
  std::cout << "Created " << lagr1_numDofs << " dofs (lagrange1).\n"
	    << "Created " << lagr0_numDofs << " dofs (lagrange0).\n"
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
  matrix_t A =  TransmissionSystemMatrix( V, K, W, alpha );


  //============================================================================
  // CREATE ANALYTICAL GRID FUNCTION FOR TRACES OF INCIDENT FIELD UINC
  //============================================================================
  typedef analytical::LinearExcitationField incident_field_t;
  typedef bem::AnalyticalGridFunction< grid_factory_t, incident_field_t,
				       Trace::Dirichlet > analytical_TDuinc_t;
  typedef bem::AnalyticalGridFunction< grid_factory_t, incident_field_t,
				       Trace::Neumann   > analytical_TNuinc_t;
  
  const incident_field_t    uinc;
  const analytical_TDuinc_t analytical_TDuinc( gridFactory, uinc );
  const analytical_TNuinc_t analytical_TNuinc( gridFactory, uinc );


  //============================================================================
  // CREATE COEFFICIENTS VECTOR FOR GRID-FUNCTIONS OF THE TRACES OF UINC
  //============================================================================
  const auto  coeff_TDui  = DofInterpolator()( analytical_TDuinc,
					       dh_lagrange1.fespace() );

  const auto  coeff_TNui  = DofInterpolator()( analytical_TNuinc,
					       dh_lagrange0.fespace() );


  //============================================================================
  // BUILD RIGHT HAND SIDE USING TRACES DEFINED ABOVE
  //============================================================================
  Eigen::VectorXd rhs(lagr0_numDofs + lagr1_numDofs,1);
  // Mass matrix with test $\qbe$ and trial $\sbe$
  typedef IdentityOperator< dh_lagrange0_t::fespace_t,
			    dh_lagrange1_t::fespace_t > id_op01_t;
  id_op01_t id_op01( dh_lagrange0.fespace(), dh_lagrange1.fespace() );
  id_op01.compute( );
  const auto& M = id_op01.matrix();
  // Evaluate <TD0 uinc, v> and <TN0 uinc, v> (taking into account that TN0 = -TN)
  rhs.segment(0            , lagr0_numDofs) = -M* coeff_TDui;
  rhs.segment(lagr0_numDofs, lagr1_numDofs) = -M.transpose()*coeff_TNui;


  //============================================================================
  // SOLVE (USING DIRECT SOLVER)
  //============================================================================
  // Define the solver
  typedef Eigen::ColPivHouseholderQR< matrix_t > directSolver_t;
  directSolver_t solverA( A );   
  // Initialize the solver
  solverA.compute( A );   
  // Solve the system
  const auto& sol = solverA.solve( rhs );   
  std::cout << "Solved linear system Ax = rhs. " << std::endl;

  std::cout << " Energy is : " << computeEnergy(sol, gridFactory, dh_lagrange0,
						dh_lagrange1)
	    << std::endl;

  return sol;

}



//------------------------------------------------------------------------------
int main( int argc, char* argv[] )
{
    Eigen::VectorXi levels(4);
    levels << 32, 128, 512, 2048;
    int Nmax = 4;
    double alpha = 2.;
    
    for(int k=0; k<Nmax; k++){
      const std::string basename ="sphere_" + std::to_string(levels(k));
      const auto sol = solveTransmissionProblem(basename, alpha);
    }

  // that's it!
  return EXIT_SUCCESS;
}


