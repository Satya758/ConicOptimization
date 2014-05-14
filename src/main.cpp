/*!
 * Tests sample data
 * FIXME Replace all these test samples and file with Google test
 */

#include <iostream>

#include <OperatorSplittingConicSolver.hpp>

using namespace std;

scs::Problem getTestData() {
  scs::Problem problem(2, 3);

  problem.b(0) = 1;
  problem.b(1) = 1;

  problem.c(0) = 1;
  problem.c(1) = 1;
  problem.c(2) = 1;

  problem.A.insert(0, 0) = 1;
  problem.A.insert(0, 1) = -1;
  problem.A.insert(1, 2) = 1;

  return problem;
}

int main(int argc, char **argv) {

  scs::Problem problem = getTestData();

  scs::ConicSolver solver(problem, 100);

  solver.solve();

  return 0;
}
