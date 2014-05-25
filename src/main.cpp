/*!
 * Tests sample data
 * FIXME Replace all these test samples and file with Google test
 */

#include <iostream>

#include <OperatorSplittingConicSolver.hpp>
#include <ProblemGeneration.hpp>

using namespace std;

scs::Problem getTestData() {
  scs::Problem problem(2, 2);

  problem.c(0) = -30;
  problem.c(1) = -40;

  problem.b(0) = 3000;
  problem.b(1) = 4000;

  problem.A.insert(0, 0) = 20;
  problem.A.insert(0, 1) = 30;

  problem.A.insert(1, 0) = 40;
  problem.A.insert(1, 1) = 30;

  return problem;
}

scs::Problem getTestData2() {
  scs::Problem problem(2, 4);

  problem.c(0) = -17.1667;
  problem.c(1) = -25.8667;
  problem.c(2) = 0;
  problem.c(3) = 0;

  problem.b(0) = 2270;
  problem.b(1) = 1900;

  problem.A.insert(0, 0) = 13;
  problem.A.insert(0, 1) = 19;
  problem.A.insert(0, 2) = 1;

  problem.A.insert(1, 0) = 20;
  problem.A.insert(1, 1) = 29;
  problem.A.insert(1, 3) = 1;

  return problem;
}

int main(int argc, char **argv) {

  scs::internal::tests::ProblemGeneration problemGeneration(400, 4, 4);

  scs::Problem problem = problemGeneration.getLinearProgram();

//   cout << problem.A << endl;
//   cout << problem.c << endl;
//   cout << problem.b << endl;

  // problem.c = problem.c * -1;
  // scs::Problem problem = getTestData2();

  scs::ConicSolver solver(problem, 2550, 1.5);

  solver.solve();


  return 0;
}
