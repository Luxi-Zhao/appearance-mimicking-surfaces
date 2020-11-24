#include "massmatrix.h"

void massmatrix(
  const Eigen::MatrixXd & l,
  const Eigen::MatrixXi & F,
  Eigen::DiagonalMatrix<double,Eigen::Dynamic> & M)
{
  // Add your code here
  int num_vertices = F.maxCoeff() + 1;
  M.resize(num_vertices);
  // Unlike SparseMatrix, resize for DiagonalMatrix doesn't
  // setZero, so call setZero separately
  M.setZero();

  for(int i = 0; i < F.rows(); i++) {
    // Get edge lengths
    double l0, l1, l2;
    l0 = l(i, 0);
    l1 = l(i, 1);
    l2 = l(i, 2);
    // Get area
    double r = 0.5 * (l0 + l1 + l2);
    double area = std::sqrt(r * (r-l0) * (r-l1) * (r-l2));
    M.diagonal()[F(i, 0)] += 1.0/3.0 * area;
    M.diagonal()[F(i, 1)] += 1.0/3.0 * area;
    M.diagonal()[F(i, 2)] += 1.0/3.0 * area;
  }
}

