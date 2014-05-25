#ifndef PRIMAL_DUAL_POINTS_HPP
#define PRIMAL_DUAL_POINTS_HPP

#include <Eigen/Dense>
#include <Eigen/Sparse>

#include "Problem.hpp"

namespace scs {
namespace internal {

// Class definition to use as friends
// TODO Too many friends is it good? Is whole point diluted? Why not make public
// if all objects are friends :-)
class Psi;
class Residuals;
class SubspaceProjection;
class ConicProjection;

// Test case
namespace tests {
class ProblemGeneration;
}

/*!
 * TODO Can we have this members variables as private?
 * Same applied to Psi object, may be create initialPoint via constructor
 */
class Omega {
 public:
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

  static Omega getFinalValue(const Omega& omega){
    Omega finalOmega;

    finalOmega.x = omega.x / omega.tau;
    finalOmega.y = omega.y / omega.tau;

    return finalOmega;
  }
  // TODO For testing only
  Eigen::VectorXd x;
 private:
  //Eigen::VectorXd x;
  Eigen::VectorXd y;
  double tau;

  friend Omega operator+(const Omega&, const Omega&);
  friend std::ostream& operator<<(std::ostream&, const Omega&);
  friend Omega operator+(const Omega&, const Psi&);
  friend Psi operator-(const Psi&, const Omega&);
  friend Omega operator-(const Omega&, const Psi&);
  friend Psi updateDualVariables(const Omega&, const Omega&, const Psi&);
  friend Omega relaxOmegaHat(const Omega&, const Omega&, const double);

  friend Psi;
  friend Residuals;
  friend SubspaceProjection;
  friend ConicProjection;
  // Test case
  friend tests::ProblemGeneration;
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

  static Psi getFinalValue(const Psi& psi, const Omega& omega){
    Psi finalPsi;

    finalPsi.s = psi.s / omega.tau;

    return finalPsi;
  }

 private:
  Eigen::VectorXd r;
  Eigen::VectorXd s;
  double kappa;

  friend std::ostream& operator<<(std::ostream&, const Psi&);
  friend Omega operator+(const Omega&, const Psi&);
  friend Psi operator-(const Psi&, const Omega&);
  friend Omega operator-(const Omega&, const Psi&);
  friend Psi updateDualVariables(const Omega&, const Omega&, const Psi&);
  friend Omega relaxOmegaHat(const Omega&, const Omega&, const double);

  friend Residuals;
  friend SubspaceProjection;
  friend ConicProjection;
  // Test case
  friend tests::ProblemGeneration;
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
 * TODO X is not relaxed...Why?
 * TODO Move outside of this class and create a as function in internal
 * namespace, but as friend to Omega or something like that
 */
Omega relaxOmegaHat(const Omega& omegaHat, const Omega& omega,
                    const double relaxationParameter) {
  Omega relaxedOmegaHat;

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
Psi updateDualVariables(const Omega& omegaHat, const Omega& omega,
                        const Psi& psi) {
  Psi updatedPsi;

  updatedPsi.r = psi.r;
  updatedPsi.s = psi.s + (omega.y - omegaHat.y);
  updatedPsi.kappa = psi.kappa + (omega.tau - omegaHat.tau);

  return updatedPsi;
}
}  // internal
}  // scs

#endif  // PRIMAL_DUAL_POINTS_HPP