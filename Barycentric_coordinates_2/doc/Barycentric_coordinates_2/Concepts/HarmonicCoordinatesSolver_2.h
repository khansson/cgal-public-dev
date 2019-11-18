/*!
\ingroup PkgBarycentric_coordinates_2Concepts
\cgalConcept

A concept that describes the set of methods that should be defined for all
solvers that can be used to solve the Laplace equation in 
`Barycentric_coordinates::Harmonic_coordinates_2`.

This concept follows the generic sparse solver concept from Eigen.
*/

class HarmonicCoordinatesSolver_2 {

public:

  /*!  
    computes factorization of a sparse square `n x n` symmetric positive definite 
    matrix `A`, where `n` is the number of unknowns.
  */
  template<typename FT>
  void compute(
    const Eigen::SparseMatrix<FT>& A) {

  }

  /*!  
    solves a linear system `Ax = b`, where `A` is a sparse square `n x n` 
    symmetric positive definite matrix, whose factorization is performed by the
    method above, `b` is a `n x 1` vector, and `n` is the number of unknowns.

    \return solution `n x 1` vector `x`.
  */
  template<typename FT>
  Eigen::Matrix<FT, Eigen::Dynamic, Eigen::Dynamic> solve(
    const Eigen::Matrix<FT, Eigen::Dynamic, Eigen::Dynamic>& b) {

  }
};
