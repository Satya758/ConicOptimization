/*!
 * TODO Name of the project is ????? not sure as of now lets see later
 * Main solver, Currently not sure how this turns out
 *
 * TODO Move code to different headers
 * TODO Change the name of the header file
 * TODO Add const to function definitions
 * @version 0.1.0
 * @author satya
 */

#ifndef OPERATOR_SPLITTING_CONIC_SOLVER_HPP
#define OPERATOR_SPLITTING_CONIC_SOLVER_HPP

#include "Logger.hpp"

// Temp
#include <iostream>

using namespace std;

#include <cmath>
#include <utility>

#include <Eigen/Sparse>
#include <Eigen/Dense>
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

// TODO How to avoid conflicting namespaces?
namespace scs {

/*!
 * TODO Currently only LP is supported, in future versions free variables, SOCP
 * constraints, SDP constraints, Exp constraints will be supported,
 * data structures for different type of constraints should be thought.
 * 1. Have all constraints in single matrix A and collect cone sizes in seperate
 *    vector.
 * 2. Collect different constraints in seperate matrices, merge them inside
 *    solver, it may be slow (I am not sure) as we have to copy them into single
 *    large matrix and collect cone sizes from the matrix sizes
 *
 * TODO Use templates to
 *
 *
 */
class Problem {
 public:
  Problem(int rows, int cols) : A(rows, cols), b(rows), c(cols) {}

  Eigen::SparseMatrix<double> A;
  Eigen::VectorXd b;
  Eigen::VectorXd c;

 private:
};

namespace internal {

/*!
 * TODO Can we have this members variables as private?
 * Same applied to Psi object, may be create initialPoint via constructor
 */
class Omega {
 public:
  Eigen::VectorXd x;
  Eigen::VectorXd y;
  double tau;

  /*!
  * TODO Move to constructor
  */
  static Omega getInitialPoint(const Problem& problem) {
    Omega initialPoint;

    initialPoint.x = Eigen::VectorXd::Zero(problem.c.rows());
    initialPoint.y = Eigen::VectorXd::Zero(problem.b.rows());
    initialPoint.tau = std::sqrt(problem.c.rows() + problem.b.rows() + 1);

    return initialPoint;
  }
};

Omega operator+(const Omega& lhsOmega, const Omega& rhsOmega) {
  Omega sumOmega;

  sumOmega.x = lhsOmega.x + rhsOmega.x;
  sumOmega.y = lhsOmega.y + rhsOmega.y;
  sumOmega.tau = lhsOmega.tau + rhsOmega.tau;

  return sumOmega;
}

std::ostream& operator<<(std::ostream& stream, const Omega& omega) {
  stream << std::endl;
  stream << "###################################" << std::endl;
  stream << "Value of x: " << std::endl << omega.x << std::endl;
  stream << "Value of y: " << std::endl << omega.y << std::endl;
  stream << "Value of tau: " << std::endl << omega.tau << std::endl;
  stream << "###################################" << std::endl;

  return stream;
}

/*!
 *
 */
class Psi {
 public:
  Eigen::VectorXd r;
  Eigen::VectorXd s;
  double kappa;

  Psi& operator=(const Psi& otherPsi) {

    r = otherPsi.r;
    s = otherPsi.s;
    kappa = otherPsi.kappa;

    return *this;
  }

  /*!
   * TODO Move to constructor
   */
  static Psi getInitialPoint(const Problem& problem) {
    Psi initialPoint;

    initialPoint.r = Eigen::VectorXd::Zero(problem.c.rows());
    initialPoint.s = Eigen::VectorXd::Zero(problem.b.rows());
    initialPoint.kappa = std::sqrt(problem.c.rows() + problem.b.rows() + 1);

    return initialPoint;
  }
};

std::ostream& operator<<(std::ostream& stream, const Psi& psi) {
  stream << std::endl;
  stream << "###################################" << std::endl;
  stream << "Value of r: " << std::endl << psi.r << std::endl;
  stream << "Value of s: " << std::endl << psi.s << std::endl;
  stream << "Value of kappa: " << std::endl << psi.kappa << std::endl;
  stream << "###################################" << std::endl;

  return stream;
}

/*!
 * FIXME If we make member variables of omega and psi private, how can we access
 * members to add them? But problem members being public is dangerous, move it
 * to private afterwards.
 *
 * Returned object is not anymore Omega but placeholder for summation, temporary
 *holder of sum instead of creating new Object.
 */

Omega operator+(const Omega& omega, const Psi& psi) {
  Omega sumOmega;

  sumOmega.x = omega.x + psi.r;
  sumOmega.y = omega.y + psi.s;
  sumOmega.tau = omega.tau + psi.kappa;

  return sumOmega;
}

/*!
 * Same description as above
 */
Psi operator-(const Psi& psi, const Omega& omega) {
  Psi minusPsi;

  minusPsi.r = psi.r - omega.x;
  minusPsi.s = psi.s - omega.y;
  minusPsi.kappa = psi.kappa - omega.tau;

  return minusPsi;
}

/*!
 * Same description as above
 */
Omega operator-(const Omega& omega, const Psi& psi) {
  Omega minusOmega;

  minusOmega.x = omega.x - psi.r;
  minusOmega.y = omega.y - psi.s;
  minusOmega.tau = omega.tau - psi.kappa;

  return minusOmega;
}

/*!
 * Calculation of Residuals
 *
 * TODO unboundedness, infeasible calculations are pending
 */
class Residuals {
 public:
  Residuals(const Problem& problem, const Omega& omega, const Psi& psi)
      : primal{getPrimal(problem, omega, psi)},
        dual{getDual(problem, omega)},
        primalDualGap{getGap(problem, omega)} {}

  Residuals(const double tolerance)
      : primal{tolerance}, dual{tolerance}, primalDualGap{tolerance} {}

  const double primal;
  const double dual;
  const double primalDualGap;

 private:
  double getPrimal(const Problem& problem, const Omega& omega, const Psi& psi) {
    double numerator =
        (problem.A * omega.x + psi.s - problem.b * omega.tau).norm();
    double denominator = omega.tau * (1 + problem.b.norm());

    return numerator / denominator;
  }

  double getDual(const Problem& problem, const Omega& omega) {
    double numerator =
        (problem.A.transpose() * omega.y + problem.c * omega.tau).norm();
    double denominator = omega.tau * (1 + problem.c.norm());

    return numerator / denominator;
  }

  double getGap(const Problem& problem, const Omega& omega) {

    double cTransposeX = problem.c.transpose() * omega.x;
    double bTranspoesY = problem.b.transpose() * omega.y;

    double numerator = std::abs(cTransposeX + bTranspoesY);
    double denominator =
        omega.tau * (1 + std::abs(cTransposeX) + std::abs(bTranspoesY));

    return numerator / denominator;
  }
};

bool operator<=(const Residuals& lhsResidual, const Residuals& rhsResidual) {
  return lhsResidual.primal <= rhsResidual.primal &&
         lhsResidual.dual <= rhsResidual.dual &&
         lhsResidual.primalDualGap <= rhsResidual.primalDualGap;
}

/*!
 * Factorization is done based on
 *http://www.princeton.edu/~rvdb/tex/myPapers/sqd6.pdf
 *
 */
class SubspaceProjection {
 public:
  SubspaceProjection(const Problem& problem, const double xScalingParameter)
      : rhsH(problem.c.rows() + problem.b.rows()), xScalingParameter(xScalingParameter) {

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

    Eigen::VectorXd deltaXTauC = delta.x * xScalingParameter - delta.tau * problem.c;
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
}  // namespace internal

/*!
*
* Interface to the solver, which accepts problem data
* TODO Currently it only supports linear programming.
*
* tolerance less than 1e-3 takes long time to converge, SCS is not meant for
* higher accuracy
*/
class ConicSolver {
 public:
  // TODO Move problem input to solve method if we dont do anything in between
  ConicSolver(const Problem& problem,
              const unsigned int maximumIterations = 1000,
              const double relaxationParameter = 1.5,
              const double xScalingParameter = 1e-3,
	      const double tolerance = 1e-3)
      : maximumIterations(maximumIterations),
        tolerance(tolerance),
        relaxationParameter(relaxationParameter),
        xScalingParameter(xScalingParameter),
        problem(problem){}

  /*!
   * TODO Add return parameter
   */
  void solve() {

    internal::Logger logger;

    int iteration = 0;

    internal::Omega omega = internal::Omega::getInitialPoint(problem);
    internal::Psi psi = internal::Psi::getInitialPoint(problem);

    // TODO How do we deal with infeasible and unboundedness
    // TODO May be we can have this in constructor
    internal::Residuals tolerantResiduals(tolerance);
    internal::SubspaceProjection subspaceProjection(problem, xScalingParameter);
    internal::ConicProjection conicProjection;

    BOOST_LOG_SEV(logger.lg, internal::logging::trivial::trace) << omega;
    while (iteration < maximumIterations) {

      internal::Residuals residuals(problem, omega, psi);
      if (residuals <= tolerantResiduals) {
        break;
      }

      // BOOST_LOG_SEV(logger.lg, internal::logging::trivial::trace) << "Initial
      // Values";
      // BOOST_LOG_SEV(logger.lg, internal::logging::trivial::trace) <<
      //"Iteration" << iteration;
      internal::Omega omegaHat =
          subspaceProjection.doProjection(problem, omega + psi);

      // BOOST_LOG_SEV(logger.lg, internal::logging::trivial::trace) <<
      // omegaHat;

      omegaHat = relaxOmegaHat(omegaHat, omega);

      // cout << "After Subspace Projection" << endl;
      // cout << omegaHat << endl;

      // BOOST_LOG_SEV(logger.lg, internal::logging::trivial::trace) <<
      // omegaHat;

      omega = conicProjection.doProjection(omegaHat - psi);

      psi = updateDualVariables(omegaHat, omega, psi);

      ++iteration;
    }
    // BOOST_LOG_SEV(logger.lg, internal::logging::trivial::trace) << "Final";

    BOOST_LOG_SEV(logger.lg, internal::logging::trivial::trace)
        << "Iteration: " << iteration;
    BOOST_LOG_SEV(logger.lg, internal::logging::trivial::trace)
        << "Final Value of X" << std::endl << omega.x / omega.tau;
  }

 private:
  const unsigned int maximumIterations;
  const double tolerance;
  const double relaxationParameter;
  const double xScalingParameter;
  const Problem& problem;


  /*!
   * TODO X is not relaxed...Why?
   * TODO Move outside of this class and create a as function in internal
   * namespace, but as friend to Omega or something like that
   */
  internal::Omega relaxOmegaHat(const internal::Omega& omegaHat,
                                const internal::Omega& omega) {
    internal::Omega relaxedOmegaHat;

    relaxedOmegaHat.x = omegaHat.x;
    relaxedOmegaHat.y =
        relaxationParameter * omegaHat.y + (1 - relaxationParameter) * omega.y;
    relaxedOmegaHat.tau = relaxationParameter * omegaHat.tau +
                          (1 - relaxationParameter) * omega.tau;

    return relaxedOmegaHat;
  }

  /*!
   * x is not relaxed/Used
   */
  internal::Psi updateDualVariables(const internal::Omega& omegaHat,
                                    const internal::Omega& omega,
                                    const internal::Psi& psi) {
    internal::Psi updatedPsi;

    updatedPsi.r = psi.r;
    updatedPsi.s = psi.s + (omega.y - omegaHat.y);
    updatedPsi.kappa = psi.kappa + (omega.tau - omegaHat.tau);

    return updatedPsi;
  }
};
}

#endif  // OPERATOR_SPLITTING_CONIC_SOLVER_HPP
