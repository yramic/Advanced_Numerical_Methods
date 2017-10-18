#ifndef INDIRECT_BEM_HPP
#define INDIRECT_BEM_HPP

#include <cmath>
#include <Eigen/Dense>

#include "../CppHilbert/Library/source/BoundaryMesh.hpp"

namespace IndirectFirstKind{

  /* 
   * @brief Build and solve indirect first kind BIE arising from Dirichlet 
   *        Laplace problem.
   * \param[in] mesh 
   * \param[in] g Dirichlet data
   * \returns coefficient vector of \f$ \mathcal{S}^{-1}_0(\mathcal{G}\f$ 
   *          corresponding to BEM solution.
   */
  template <typename FUNC>
  Eigen::VectorXd solve(const BoundaryMesh& mesh, const FUNC& g){


  }


}



namespace IndirectSecondKind{

  /* 
   * @brief Build and solve indirect second kind BIE arising from Dirichlet 
   *        Laplace problem. 
   * \param[in] mesh 
   * \param[in] g Dirichlet data
   * \returns coefficient vector of \f$ \mathcal{S}^{-1}_0(\mathcal{G}\f$ 
   *          corresponding to BEM solution.
   */
  template <typename FUNC>
  Eigen::VectorXd solve(const BoundaryMesh& mesh, const FUNC& g){


  }


}


#endif
