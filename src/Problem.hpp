#ifndef PROBLEM_HPP
#define PROBLEM_HPP

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
 * TODO Use templates to ??? what??
 * TODO Add configuration information (log info, solver control parameters) and
 *output structures like info,
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
}

#endif  // PROBLEM_HPP