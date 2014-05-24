#ifndef PROJECTIONS_HPP
#define PROJECTIONS_HPP

#include <Eigen/Dense>
#include <Eigen/Sparse>

//****************************************************************************//
// TODO Why do I have to add this again they were any way added in
// <Eigen/PaStiXSupport>, but are not considered for some reason, without these
// there are compiler errors
#include <complex>
extern "C" {
#include <pastix_nompi.h>
#include <pastix.h>
}

#define COMPLEX std::complex<float>
#define DCOMPLEX std::complex<double>
//****************************************************************************//

#include <Eigen/PaStiXSupport>

#include "Problem.hpp"
#include "PrimalDualPoints.hpp"

namespace scs {
namespace internal {

/*!
* Factorization is done based on
*http://www.princeton.edu/~rvdb/tex/myPapers/sqd6.pdf
*
*/
class SubspaceProjection {
 public:
  SubspaceProjection(const Problem& problem, const double xScalingParameter)
      : rhsH(problem.c.rows() + problem.b.rows()),
        xScalingParameter(xScalingParameter) {

    Eigen::SparseMatrix<double> MHat = createMHat(problem);
    // Factorize
    directSolver.compute(MHat);
    if (directSolver.info() != Eigen::Success) {
      // TODO Exception or what? BOOST Log?
      // TODO How to have global logger object...singleton
      cout << "Failed" << endl;
    }

    rhsH << problem.c, problem.b;

    MInverseH = directSolver.solve(rhsH);
    // Negate y part as we have normalized M to MHat to confirm with spd
    MInverseH.tail(problem.b.rows()) = -MInverseH.tail(problem.b.rows());

    // Denominator of Matrix inverse lemma
    denominator = 1 + rhsH.transpose() * MInverseH;
  }

  /*!
   *
   */
  Omega doProjection(const Problem& problem, const Omega& delta) const {

    Eigen::VectorXd deltaXTauC =
        delta.x * xScalingParameter - delta.tau * problem.c;
    Eigen::VectorXd deltaYTauB = delta.y - delta.tau * problem.b;
    Eigen::VectorXd deltaHat(problem.c.rows() + problem.b.rows());
    deltaHat << deltaXTauC, deltaYTauB;

    Eigen::VectorXd MInverseDelta = directSolver.solve(deltaHat);
    // Negate y part as we have normalized M to MHat to confirm with spd
    MInverseDelta.tail(problem.b.rows()) =
        -MInverseDelta.tail(problem.b.rows());

    // (mil) Matrix Inverse Lemma solution
    Eigen::VectorXd milSolution =
        MInverseDelta -
        (MInverseH * rhsH.transpose() * MInverseDelta) / denominator;

    Omega omegaHat;
    omegaHat.x = milSolution.head(problem.A.cols());
    omegaHat.y = milSolution.tail(problem.A.rows());
    omegaHat.tau = delta.tau + problem.c.transpose() * omegaHat.x +
                   problem.b.transpose() * omegaHat.y;

    return omegaHat;
  }

 private:
  Eigen::PastixLDLT<Eigen::SparseMatrix<double>, Eigen::Lower> directSolver;
  // Both First part and second part
  Eigen::VectorXd MInverseH;
  Eigen::VectorXd rhsH;
  // Denominator for Matrix Inverse lemma
  double denominator;

  const double xScalingParameter;

  /*!
   * Create MHat which confirms to SPD paper
   * TODO Is there a better way to do rather than working with column compressed
   * format
   */
  Eigen::SparseMatrix<double> createMHat(const Problem& problem) {

    int Arows = problem.A.rows();
    int Acols = problem.A.cols();
    int squareMatrixSize = Arows + Acols;

    Eigen::SparseMatrix<double> lowerMHat(squareMatrixSize, squareMatrixSize);

    // As storage schem is column compressed, we fill data column wise to
    // achieve O(1) complexity. And we know it is symmetric so we fill only
    // lower traingular matrix
    for (int colIndex = 0; colIndex < squareMatrixSize; ++colIndex) {
      // Change the rowIndex as we need to traverse only lower triangle, thats
      // the reason for rowIndex = colIndex
      for (int rowIndex = colIndex; rowIndex < squareMatrixSize; ++rowIndex) {
        // fill top left blocks diagonal
        if (colIndex == rowIndex && colIndex < Acols && rowIndex < Acols) {
          lowerMHat.insert(rowIndex, colIndex) = xScalingParameter;
        }
        // fill bottom right block diagonal
        if (colIndex == rowIndex && colIndex >= Acols && rowIndex >= Acols) {
          lowerMHat.insert(rowIndex, colIndex) = -1;
        }
        // fill bottom left block using A
        if (rowIndex >= Acols && colIndex < Acols) {
          lowerMHat.insert(rowIndex, colIndex) =
              -problem.A.coeff(rowIndex - Acols, colIndex);
        }
      }
    }

    return lowerMHat;
  }
};

/*!
 *
 */
class ConicProjection {
 public:
  /*!
   * TODO As of now its only LP cone
   */
  Omega doProjection(const Omega& omegaHat) {

    Omega omegaHat2 = omegaHat;
    omegaHat2.y = omegaHat2.y.unaryExpr([](const double & element)->double {
      if (element <= 0) {
        return 0;
      } else {
        return element;
      }
    });

    if (omegaHat2.tau < 0) {
      omegaHat2.tau = 0;
    }

    return omegaHat2;
  }
};

}  // internal
}  // scs
#endif  // PROJECTIONS_HPP