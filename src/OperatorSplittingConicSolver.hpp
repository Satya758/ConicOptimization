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

// Temp
#include <iostream>

using namespace std;

#include <Eigen/Sparse>
#include <Eigen/Dense>

#include "Problem.hpp"

#include "Core/Logger.hpp"
#include "Core/PrimalDualPoints.hpp"
#include "Core/Residuals.hpp"
#include "Core/Projections.hpp"

// TODO How to avoid conflicting namespaces?
namespace scs {

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
              const double tolerance = 1e-6)
      : maximumIterations(maximumIterations),
        tolerance(tolerance),
        relaxationParameter(relaxationParameter),
        xScalingParameter(xScalingParameter),
        problem(problem) {}

  /*!
   * TODO Add return parameter
   */
  void solve() {

    internal::Logger logger;

    internal::Omega omega = internal::Omega::getInitialPoint(problem);
    internal::Psi psi = internal::Psi::getInitialPoint(problem);

    // TODO How do we deal with infeasible and unboundedness
    // TODO May be we can have this in constructor
    internal::Residuals tolerantResiduals(tolerance);
    internal::SubspaceProjection subspaceProjection(problem, xScalingParameter);
    internal::ConicProjection conicProjection;

    SolverState solverState;
    // TODO Move it to for loop and information to info object which is return parameter of solver method
    int iteration;
    for (iteration = 0; iteration < maximumIterations; ++iteration) {

      internal::Residuals residuals(problem, omega, psi);

      solverState = internal::getSolverState(residuals, tolerantResiduals);

      if (solverState == SolverState::CONVERGED ||
          solverState == SolverState::UNBOUNDED ||
          solverState == SolverState::INFEASIBLE) {
        break;
      }

      internal::Omega omegaHat =
          subspaceProjection.doProjection(problem, omega + psi);

      omegaHat = internal::relaxOmegaHat(omegaHat, omega, relaxationParameter);

      omega = conicProjection.doProjection(omegaHat - psi);

      psi = internal::updateDualVariables(omegaHat, omega, psi);
    }

    cout << "SolverState: " << internal::getSolverState(solverState) << endl;
        BOOST_LOG_SEV(logger.lg, internal::logging::trivial::trace)
            << "Iteration: " << iteration;
    /*BOOST_LOG_SEV(logger.lg, internal::logging::trivial::trace)
        << "Final Value of X" << std::endl << internal::Omega::getFinalValue(omega)*/;
	BOOST_LOG_SEV(logger.lg, internal::logging::trivial::trace)
            << "Solver Objective value: " << problem.c.transpose() * internal::Omega::getFinalValue(omega).x;

  }

 private:
  const unsigned int maximumIterations;
  const double tolerance;
  const double relaxationParameter;
  const double xScalingParameter;
  const Problem& problem;
};
}

#endif  // OPERATOR_SPLITTING_CONIC_SOLVER_HPP
