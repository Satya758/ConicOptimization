#ifndef RESIDUALS_HPP
#define RESIDUALS_HPP

#include <cmath>
#include <string>

#include <Eigen/Dense>
#include <Eigen/Sparse>

#include "Problem.hpp"

#include "PrimalDualPoints.hpp"

namespace scs {
namespace internal {

/*!
 * Calculation of Residuals
 *
 * TODO unboundedness, infeasible calculations are pending
 */
class Residuals {
 private:
  double normOfC;
  double normOfB;
  Eigen::VectorXd AXPlusS;
  Eigen::VectorXd ATransposeY;
  double cTransposeX;
  double bTranspoesY;

  double getPrimal(const Problem& problem, const Omega& omega,
                   const Psi& psi) const {
    double numerator = (AXPlusS - problem.b * omega.tau).norm();
    double denominator = omega.tau * (1 + normOfB);

    return numerator / denominator;
  }

  double getDual(const Problem& problem, const Omega& omega) const {
    double numerator = (ATransposeY + problem.c * omega.tau).norm();
    double denominator = omega.tau * (1 + normOfC);

    return numerator / denominator;
  }

  double getGap(const Problem& problem, const Omega& omega) const {

    double numerator = std::abs(cTransposeX + bTranspoesY);
    double denominator =
        omega.tau + std::abs(cTransposeX) + std::abs(bTranspoesY);

    return numerator / denominator;
  }

 public:
  Residuals(const Problem& problem, const Omega& omega, const Psi& psi)
      : normOfC{problem.c.norm()},
        normOfB{problem.b.norm()},
        AXPlusS{problem.A * omega.x + psi.s},
        ATransposeY{problem.A.transpose() * omega.y},
        cTransposeX{problem.c.transpose() * omega.x},
        bTranspoesY{problem.b.transpose() * omega.y},
        primal{getPrimal(problem, omega, psi)},
        dual{getDual(problem, omega)},
        primalDualGap{getGap(problem, omega)},
        unbounded{cTransposeX < 0 ? AXPlusS.norm() * normOfC / -cTransposeX
                                  : NAN},
        infeasible{bTranspoesY < 0 ? ATransposeY.norm() * normOfB / -bTranspoesY
                                   : NAN} {}

  Residuals(const double tolerance)
      : primal{tolerance},
        dual{tolerance},
        primalDualGap{tolerance},
        unbounded{tolerance},
        infeasible{tolerance} {}

  const double primal;
  const double dual;
  const double primalDualGap;
  const double unbounded;
  const double infeasible;
};

SolverState getSolverState(const Residuals& currentResidual,
                           const Residuals& tolerance) {

  if (currentResidual.unbounded <= tolerance.unbounded) {
    return SolverState::UNBOUNDED;
  } else if (currentResidual.infeasible <= tolerance.infeasible) {
    return SolverState::INFEASIBLE;
  } else if (currentResidual.primal < tolerance.primal &&
             currentResidual.dual < tolerance.dual &&
             currentResidual.primalDualGap < tolerance.primalDualGap) {
    return SolverState::CONVERGED;
  } else {
    return SolverState::IN_PROGRESS;
  }
}

/*!
 * TODO Used in Logging, May be move it to logger class!?
 *
 */
std::string getSolverState(const SolverState solverState) {

  switch (solverState) {
    case SolverState::CONVERGED:
      return "Converged";
    case SolverState::INFEASIBLE:
      return "Infeasible";
    case SolverState::UNBOUNDED:
      return "Unbounded";
    case SolverState::IN_PROGRESS:
      return "In Progress";
  }
}

}  // internal
}  // scs

#endif  // RESIDUALS_HPP