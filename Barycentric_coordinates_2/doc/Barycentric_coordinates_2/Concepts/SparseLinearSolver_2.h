namespace CGAL {
namespace Barycentric_coordinates {

/*!
\ingroup PkgBarycentricCoordinates2RefConcepts
\cgalConcept

A concept that describes the set of methods that should be defined for all
sparse linear solvers that can be used to solve the Laplace equation.

This concept follows the generic sparse solver concept from \ref thirdpartyEigen "Eigen".

\cgalHasModel All models of \ref thirdpartyEigen "Eigen" sparse solvers.
*/

class SparseLinearSolver_2 {

public:

    /// \name Types
    /// @{
      
	  /// A model of `FieldNumberType`. 
    typedef unspecified_type FT;

    /// An \ref thirdpartyEigen "Eigen" matrix type.
    typedef Eigen::SparseMatrix<FT> MatrixFT;

    /// An \ref thirdpartyEigen "Eigen" vector type.
    typedef Eigen::Matrix<FT, Eigen::Dynamic, Eigen::Dynamic> VectorFT;

    /// @}

  /*!  
    computes factorization of a sparse square `n x n` symmetric positive definite 
    matrix `A`, where `n` is the number of unknowns.
  */
  void compute(
    const MatrixFT& A) {

  }

  /*!  
    solves a linear system `Ax = b`, where `A` is a sparse square `n x n` 
    symmetric positive definite matrix, whose factorization is performed by the
    method above, `b` is a `n x 1` vector, and `n` is the number of unknowns.

    \return solution `n x 1` vector `x`.
  */
  VectorFT solve(
    const VectorFT& b) {

  }
};
}
}
