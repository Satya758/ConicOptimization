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

#include <cmath>
#include <utility>

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

      omegaHat = relaxOmegaHat(omegaHat, omega, relaxationParameter);

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
        << "Final Value of X" << std::endl << omega;
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
