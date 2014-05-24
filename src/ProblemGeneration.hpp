/*!
 *
 * Move to Google Test folder/project later
 *
 */
#ifndef PROBLEM_GENERATION_HPP
#define PROBLEM_GENERATION_HPP

#include <Eigen/Dense>
#include <Eigen/Sparse>

#include <OperatorSplittingConicSolver.hpp>

namespace scs {
namespace internal {
namespace tests {

// Temp TODO
#include <iostream>
using namespace std;

/*!
 * TODO Return problem object
 */
scs::Problem getLPProblem(int rows, int cols) {

  Eigen::VectorXd x(cols);
  x.setRandom(cols);

  Eigen::MatrixXd ADense = Eigen::MatrixXd::Random(rows, cols).cwiseAbs();

  scs::Problem problem(rows, cols);
  problem.A = ADense.sparseView();

  problem.c.setRandom(cols);

  problem.b = problem.A * x;

  return problem;
}
}  // tests
}  // internal
}  // scs

#endif  // PROBLEM_GENERATION_HPP