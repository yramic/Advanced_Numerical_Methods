//// 
//// Copyright (C) 2017 SAM (D-MATH) @ ETH Zurich
//// Author(s): curzuato < > 
//// Contributors:  dcasati 
//// This file is part of the AdvNumCSE repository.
////
/***********************************************************************
 *                                                                     *
 * Code for Course "Advanced Numerical Methods for CSE"                *
 * (Prof. Dr. R. Hiptmair)                                             * 
 * Author:                                                             *
 * Date:                                                               *
 * (C) Seminar for Applied Mathematics, ETH Zurich                     *
 * This code can be freely used for non-commercial purposes as long    *
 * as this header is left intact.                                      *
 ***********************************************************************/

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
// bem operator ----------------------------------------------------------------
#include <bem_operator/bem_operator.hpp>
// sparse operators ------------------------------------------------------------
#include <sparse_operators/identity_operator.hpp>
#include <sparse_operators/combinatorial_gradient.hpp>


using namespace betl2;
namespace big = betl2::input::gmsh;


//------------------------------------------------------------------------------
int main( int argc, char* argv[] )
{
  //============================================================================
  // READ MESH
  //============================================================================
  // Parse the command line to get name of mesh
  const std::string basename = betl2::parseCommandLine( argc, argv );
  
  // Create input from it
  std::cout << "Input from: " << basename << ".msh" << std::endl;
  big::Input input( basename );
  
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
  // Defining the discrete FE-spaces requires: Choosing the basis and creating 
  // the dof hanfler for it ( see NPDE ).
  // 1 - Define order of approximation. Here we choose Linear (constant and 
  // quadratic are also available).
  const fe::ApproxOrder order = fe::Linear;
  // Define bases for $\qbe, \sbe$ and edge functions (to be used when implementing
  // the Hypersingular operator using integration by parts formula (1.3.92))
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
  // Define the fundamental solution to build the BEM operators. 
  typedef bem::FundSol< bem::FSType::LAPLACE, 3 > laplace_fs_t;
  laplace_fs_t laplace_fs;

  
  //============================================================================
  // GALERKIN KERNELS
  //============================================================================
  // Define the kernel of the BEM operator using the fundamental solution and 
  // the spaces to be used. We also need to specify the underlaying layer
  // potential:
  // - SL: Single Layer,  or
  // - DL : Double Layer

  // We need to define 3 GalerkinKernel-objects to build the 4 BEM operators:
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
  // Define the singularity detector in order to identify the different singularity
  // situations on the given grid (as discussed in section 1.5.3: coinciding panels,
  // adjacent panels, common vertex and panels with positive distance). 
  typedef bem::GalerkinSingularityDetector<grid_factory_t> singularity_detector_t;
  singularity_detector_t singularity_detector( gridFactory );


  //============================================================================
  // DEFINE QUADRATURE
  //============================================================================
  // Define the quadrature rule by setting number of points for each case (positive
  // distance and singular cases)
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
  // The following boolean value is optional. If not given it is assumed to be
  // true. It indicates whether the gram determinants are to be included into the 
  // integral or not.
  // In general, these determinants are mandatory. However, certain linear
  // combinations of scalar potentials might be used to express potentials related
  // to edge functions. In that particular case, the gram determinants cancel out
  // and should(must!) be neglected (but this is not the case here, so set it true).
  const bool withGram = true;
  // Define integrator to be used with the kernels above.
  typedef bem::IntegrationTraits< quadrature_rule_t,
                                  grid_factory_t, 
                                  singularity_detector_t, 
                                  cache_t, 
                                  withGram, /* on element X */
                                  withGram  /* on element Y */ > integrationTraits;

  typedef bem::GalerkinIntegrator< lagr_sl_kernel_t,
				   integrationTraits > lagr_sl_integrator_t;
  typedef bem::GalerkinIntegrator< lagr_dl_kernel_t,
				   integrationTraits > lagr_dl_integrator_t;
  // An integrator for edge functions: In this case the withGram-flag within
  // 'integrationTraits' has no effect since integrators for edge-functions neglect
  // the Gram determinants in any case.
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
  // We create our BEM operators using the integrators defined above and the
  // required FE-spaces. We retrieve these spaces from the dofhandler constructed
  // for each one of them.
  
  // - V 
  typedef BemOperator< lagr_sl_integrator_t,
                       typename dh_lagrange0_t::fespace_t > bem_op_V_t;
  bem_op_V_t   bem_op_V  ( lagr_sl_integrator , dh_lagrange0.fespace() );
  sec_timer_t timer_op_V;
  bem_op_V.compute( );
  std::cout << "Computation of <V phi, phi>, phi in H^{-1/2} took: "
	    << timer_op_V << std::endl;
  const auto& V    = bem_op_V.matrix();
  std::cout<< "Dim check: V is " << V.rows() << " x " << V.cols() << std::endl;

  // - K
  typedef BemOperator< lagr_dl_integrator_t,
                       typename dh_lagrange0_t::fespace_t,
                       typename dh_lagrange1_t::fespace_t > bem_op_K_t;
  bem_op_K_t   bem_op_K  ( lagr_dl_integrator , dh_lagrange0.fespace(),
			   dh_lagrange1.fespace() );
  sec_timer_t timer_op_K;
  bem_op_K.compute( );
  std::cout << "Computation of <K psi, phi>, phi in H^{-1/2}, psi in H^{1/2} took: " 
            << timer_op_K << std::endl;
  const auto& K    = bem_op_K.matrix();
  std::cout<< "Dim check: K is " << K.rows() << " x " << K.cols() << std::endl;

  // - W (using integration by parts):
  // (i) Single layer operator with div basis functions
  typedef BemOperator< div_sl_integrator_t,
                       typename dh_div_t::fespace_t > bem_op_div_t;
  bem_op_div_t bem_op_div( div_sl_integrator, dh_div.fespace() );
  sec_timer_t timer_op_div;
  bem_op_div.compute( );
  std::cout << "Computation of <V u, u>, u in H^{-1/2}(div) took: "
	    << timer_op_div << std::endl;
  const auto& Vdiv = bem_op_div.matrix();

  // (ii) Discrete embedding from H(lagrange) to H(div) ( representing grad X n 
  // on formula 1.3.92).
  typedef CombinatorialGradient< typename dh_lagrange1_t::fespace_t,
                                 typename dh_div_t      ::fespace_t > cgrad_op_t;
  cgrad_op_t cgrad_op( dh_lagrange1.fespace(), dh_div.fespace() );
  cgrad_op.compute( cache );
  const auto& C = cgrad_op.matrix();

  // (iii) Then we combine these ingredients to build the matrix for W:
  const auto& W = C * Vdiv * C.transpose();
  std::cout<< "Dim check: W is " << W.rows() << " x " << W.cols() << std::endl;

  
  //============================================================================
  // CREATE MASS-MATRIX
  //============================================================================
  // We can create the mass matrix $\qbe/ \sbe$ via the IdentityOperator, by
  // specifying these FE-spaces
  typedef IdentityOperator< typename dh_lagrange0_t::fespace_t,
                            typename dh_lagrange1_t::fespace_t > discrete_op_t;
  discrete_op_t discrete_op( dh_lagrange0.fespace(), dh_lagrange1.fespace() );
  discrete_op.compute( );
  const auto& M = discrete_op.matrix();
  std::cout<< "Dim check: M is " << M.rows() << " x " << M.cols() << std::endl;


  // that's it!
  std::cout << "Done! " << std::endl;
  return EXIT_SUCCESS;
}


