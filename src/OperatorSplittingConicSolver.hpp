/*!
 * TODO Name of the project is ????? not sure as of now lets see later
 * Main solver, Currently not sure how this turns out
 * @version 0.1.0
 * @author satya
 */

// TODO For testing
#include <iostream>
using std::cout;
using std::endl;

#include <cmath>
#include <ratio>

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

#ifndef OPERATOR_SPLITTING_CONIC_SOLVER_HPP
#define OPERATOR_SPLITTING_CONIC_SOLVER_HPP

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
  cout << "Primal: " << lhsResidual.primal << endl;
  cout << "Dual: " << lhsResidual.dual << endl;
  cout << "Gap: " << lhsResidual.primalDualGap << endl;

  return lhsResidual.primal <= rhsResidual.primal &&
         lhsResidual.dual <= rhsResidual.dual &&
         lhsResidual.primalDualGap <= rhsResidual.primalDualGap;
}

/*!
 *
 */
class SubspaceProjection {
 public:
  SubspaceProjection(const Problem& problem) {
    // Factorize
    directSolver.compute(getNormalEquationLhs(problem));
    if (directSolver.info() != Eigen::Success) {
      // TODO Exception or what? BOOST Log?
      cout << "Factorization failed!" << endl;
    }

    // Solve
    // Instead of merging h, we have seperated h into MInverseHC & MInverseHB
    Eigen::VectorXd rhsH = problem.c - problem.A.transpose() * problem.b;
    MInverseHC = directSolver.solve(rhsH);
    MInverseHB = problem.b + problem.A * MInverseHC;

    // Denominator of Matrix inverse lemma
    denominator = 1 + problem.c.transpose() * MInverseHC +
                  problem.b.transpose() * MInverseHB;
  }

  /*!
   *
   */
  Omega doProjection(const Problem& problem, const Omega& delta) const {

    Eigen::VectorXd deltaXTauC = delta.x - delta.tau * problem.c;
    Eigen::VectorXd deltaYTauB = delta.y - delta.tau * problem.b;

    Eigen::VectorXd MInverseDeltaX =
        directSolver.solve(deltaXTauC - problem.A.transpose() * deltaYTauB);
    Eigen::VectorXd MInverseDeltaY = deltaYTauB + problem.A * MInverseDeltaX;

    // For some reason CT*x + bT*Y expression is failing, so the following
    // workaround
    double CTransposeDeltaX = problem.c.transpose() * MInverseDeltaX;
    double hTransposeMInverseDelta =
        CTransposeDeltaX + problem.b.transpose() * MInverseDeltaY;

    Omega omegaHat;
    omegaHat.x =
        MInverseDeltaX - MInverseHC * (hTransposeMInverseDelta / denominator);
    omegaHat.y =
        MInverseDeltaY - MInverseHB * (hTransposeMInverseDelta / denominator);
    omegaHat.tau = delta.tau + problem.c.transpose() * omegaHat.x +
                   problem.b.transpose() * omegaHat.y;

    return omegaHat;
  }

 private:
  Eigen::PastixLLT<Eigen::SparseMatrix<double>, Eigen::Lower> directSolver;
  // Both First part and second part
  Eigen::VectorXd MInverseHC;
  Eigen::VectorXd MInverseHB;
  // Denominator for Matrix Inverse lemma
  double denominator;

  Eigen::SparseMatrix<double> getNormalEquationLhs(const Problem& problem) {
    Eigen::SparseMatrix<double> identity(problem.A.cols(), problem.A.cols());
    identity.setIdentity();

    Eigen::SparseMatrix<double> ATransposeA = problem.A.transpose() * problem.A;

    return identity + ATransposeA;
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
    omegaHat2.y.unaryExpr([](const double & element)->double {
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
  ConicSolver(Problem& problem, const unsigned int maximumIterations = 1000,
              const double tolerance = 1e-3,
              const double relaxationParameter = 1.5)
      : maximumIterations(maximumIterations),
        tolerance(tolerance),
        relaxationParameter(relaxationParameter),
        problem(problem) {}

  void solve() {
    int iteration = 0;

    internal::Omega omega = internal::Omega::getInitialPoint(problem);
    internal::Psi psi = internal::Psi::getInitialPoint(problem);

    // TODO How do we deal with infeasible and unboundedness
    // TODO May be we can have this in constructor
    internal::Residuals tolerantResiduals(tolerance);
    internal::SubspaceProjection subspaceProjection(problem);
    internal::ConicProjection conicProjection;

    internal::Omega omegaHat;

    while (iteration < maximumIterations) {

      cout << "Final Result: " << omega.x << endl;
      cout << "Final Result: " << omega.y << endl;
      cout << "Final Result: " << omega.tau << endl;

      internal::Residuals residuals(problem, omega, psi);
      if (residuals <= tolerantResiduals) {
        break;
      }
      omegaHat = subspaceProjection.doProjection(problem, omega);
      omega = conicProjection.doProjection(omegaHat - psi);
      psi = psi - (omegaHat + omega);

      // cout << "Gap: " << residuals.primalDualGap << endl;

      ++iteration;
    }


  }

 private:
  const unsigned int maximumIterations;
  const double tolerance;
  const double relaxationParameter;

  const Problem& problem;
};
}

#endif  // OPERATOR_SPLITTING_CONIC_SOLVER_HPP
