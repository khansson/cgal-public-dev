/*!
\ingroup PkgBarycentric_coordinates_2Concepts
\cgalConcept

A concept that describes the set of methods that should be defined for all
solvers that can be used to solve the Laplace equation for harmonic coordinates.
*/

class HarmonicCoordinatesSolver {

public:

  /*!  
    computes factorization of the matrix A.
  */
  void compute(const Eigen::SparseMatrix<double>& A) {

  }

  /*!  
    solves a linear system `Ax = b`.

    \return solution matrix `x`.
  */
  Eigen::MatrixXd solve(const Eigen::MatrixXd& b) {

  }
};
