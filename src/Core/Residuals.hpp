#ifndef RESIDUALS_HPP
#define RESIDUALS_HPP

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
}  // internal
}  // scs

#endif  // RESIDUALS_HPP