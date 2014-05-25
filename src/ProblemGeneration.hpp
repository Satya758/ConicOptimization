/*!
 *
 * Move to Google Test folder/project later
 *
 */
#ifndef PROBLEM_GENERATION_HPP
#define PROBLEM_GENERATION_HPP

#include <random>

#include <Eigen/Dense>
#include <Eigen/Sparse>

#include "Problem.hpp"
#include "Core/PrimalDualPoints.hpp"
#include "Core/Projections.hpp"

namespace scs {
namespace internal {
namespace tests {

// Temp TODO
#include <iostream>
using namespace std;

class ProblemGeneration {
 public:
  ProblemGeneration(const int rows, const int cols, const int colNNZ)
      : rows{rows},
        cols{cols},
        colNNZ{colNNZ},
        generator{randomDevice()},
        rowDistribution{0, rows - 1},
        coeffDistribution{-10, 10} {}

  /*!
  *
  */
  scs::Problem getLinearProgram() {

    scs::Problem problem(rows, cols);

    problem.A = generateSparseMatrix();

    // TODO Use Omega and Psi instead of individual vectors, but Omega & Psi
    // constructors have to changed to accept Problem as input so that they can
    // create vector of problem size, this change is applicable to solver not
    // only to test case
    Eigen::VectorXd z(rows);

    scs::internal::Omega omega;
    scs::internal::Psi psi;

    omega.x = generateVector(cols);
    omega.y = z = generateVector(rows);

    scs::internal::ConicProjection conicProjection;
    omega = conicProjection.doProjection(omega);

    psi.s = omega.y - z;

    problem.c = -problem.A.transpose() * omega.y;
    problem.b = problem.A * omega.x + psi.s;

    cout << "Objective value: " << problem.c.transpose() * omega.x << endl;

    return problem;
  }

 private:
  const int rows;
  const int cols;
  const int colNNZ;
  std::random_device randomDevice;
  std::mt19937_64 generator;
  std::uniform_int_distribution<int> rowDistribution;
  std::uniform_real_distribution<double> coeffDistribution;

  /*!
   *
   * colNNZ is number non zero elements in column
   */
  Eigen::SparseMatrix<double> generateSparseMatrix() {

    Eigen::SparseMatrix<double> A(rows, cols);

    const int NNZ = 2 * cols * colNNZ;

    A.reserve(NNZ);

    for (int col = 0; col < cols; ++col) {
      for (int row = 0; row < colNNZ; ++row) {
        A.coeffRef(rowDistribution(generator), col) =
            coeffDistribution(generator);
      }
    }

    return A;
  }

  Eigen::VectorXd generateVector(const int size) {

    Eigen::VectorXd vector(size);

    vector = vector.unaryExpr([&](const double & element)->double {
      return coeffDistribution(generator);
    });

    return vector;
  }
};

}  // tests
}  // internal
}  // scs

#endif  // PROBLEM_GENERATION_HPP