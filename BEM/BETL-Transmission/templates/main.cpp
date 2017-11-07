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
#include <analytical_functions/linear_excitation_field.hpp>


using namespace betl2;
namespace big = betl2::input::gmsh;

typedef Eigen::MatrixXd matrix_t;
const std::string path = "../BEM/BETL-Transmission/meshes/";


//------------------------------------------------------------------------------
/* SAM_LISTING_BEGIN_0 */
matrix_t TransmissionSystemMatrix( const matrix_t&  V, const matrix_t&  K,
				   const matrix_t&  W, const double & alpha){
  matrix_t A(K.rows() + W.rows(), K.cols() + V.cols());
  // TODO: Compute block matrix A
  return A;
}
/* SAM_LISTING_END_0 */


//------------------------------------------------------------------------------
/* SAM_LISTING_BEGIN_2 */
template<typename GRID_FACTORY, typename DH_LAGR0, typename DH_LAGR1>
double computeEnergy(const Eigen::VectorXd& sol, const GRID_FACTORY gridFactory,
		     const DH_LAGR0& dh_lagr0, const DH_LAGR1& dh_lagr1){
  // TODO: Implement your code
  return 0.;
}
/* SAM_LISTING_END_2 */


//------------------------------------------------------------------------------
/* SAM_LISTING_BEGIN_1 */
Eigen::VectorXd solveTransmissionProblem(const std::string meshname,
					 const double& alpha){
  // TODO: Implement your code
 Eigen:;VectorXd sol;
  return sol;
}
/* SAM_LISTING_END_1 */



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


