#include "smooth.h"
#include <igl/edge_lengths.h>
#include <igl/cotmatrix.h>
#include <igl/massmatrix.h>
#include <igl/repdiag.h>

typedef Eigen::Triplet<double> T;

// D_A
void voronoi_area(
  const Eigen::MatrixXd & V,
  const Eigen::MatrixXi & F,
  Eigen::DiagonalMatrix<double, Eigen::Dynamic> & D_A)
{
  Eigen::SparseMatrix<double> M;
  igl::massmatrix(V, F, igl::MASSMATRIX_TYPE_VORONOI, M);
  Eigen::VectorXd M_diag = M.diagonal();
  assert(M_diag.size() == V.rows());

  D_A.resize(V.rows() * 3);
  D_A.setZero();

  for(int i = 0; i < M_diag.size(); i++) {
    for(int j = 0; j < 3; j++) {
      D_A.diagonal()[3 * i + j] = sqrt(M_diag[i]);
    }
  }
}

// L_tilda
void laplacian(
  const Eigen::MatrixXd & V,
  const Eigen::MatrixXi & F,
  Eigen::SparseMatrix<double> & L_tilda)
{
  Eigen::SparseMatrix<double> L;
  igl::cotmatrix(V, F, L);

  std::vector<T> tripletList;
  tripletList.reserve(V.rows() * V.rows() * 3);

  for(int i = 0; i < V.rows(); i++) {
    for(int j = 0; j < V.rows(); j++) {
      double l = L.coeff(i, j);
      for(int k = 0; k < 3; k++) {
        tripletList.push_back(T(3 * i + k, 3 * j + k, l));
      }
    }
  }
  L_tilda.resize(3 * V.rows(), 3 * V.rows());
  L_tilda.setFromTriplets(tripletList.begin(), tripletList.end());
}

// D_v
void sparse_v(
  const Eigen::MatrixXd & V,
  Eigen::DiagonalMatrix<double, Eigen::Dynamic> & D_v)
{
  D_v.resize(V.rows() * 3);
  D_v.setZero();

  for(int i = 0; i < V.rows(); i++) {
    for(int j = 0; j < 3; j++) {
      D_v.diagonal()[3 * i + j] = V(i, j);
    }
  }
}

// S
void selector(
  const Eigen::MatrixXd & V,
  Eigen::SparseMatrix<double> & S)
{
  Eigen::SparseMatrix<double> ones(3, 0);
  for(int i = 0; i < 3; i++) {
    ones.insert(i, 0) = 1.0;
  }

  S.resize(3 * V.rows(), V.rows());
  igl::repdiag(ones, V.rows(), S);
}

void get_lambda(
  const Eigen::MatrixXd & V,
  const Eigen::RowVector3d & o,
  Eigen::VectorXd & lambda)
{
  lambda.resize(V.rows());
  lambda << (V.rowwise() - o).norm();
}

void get_L_theta(
  const Eigen::VectorXd & S_lambda0,
  const Eigen::SparseMatrix<double> L_tilda0,
  const Eigen::DiagonalMatrix<double, Eigen::Dynamic> D_v,
  Eigen::VectorXd & L_theta)
{
  Eigen::DiagonalMatrix<double, Eigen::Dynamic> D_S_lambda0;
  D_S_lambda0.diagonal() = S_lambda0;

  L_theta.setZero();
  L_theta = D_S_lambda0.inverse() * L_tilda0 * D_v * S_lambda0;
}

void deform(
  const Eigen::MatrixXd & V,
  const Eigen::MatrixXi & F,
  Eigen::VectorXd & lambda)
{
  // D_A
  Eigen::DiagonalMatrix<double, Eigen::Dynamic> D_A;
  voronoi_area(V, F, D_A);

  // S
  Eigen::SparseMatrix<double> S;
  selector(V, S);

  // lambda0
  Eigen::RowVector3d o;
  o.setZero();
  Eigen::VectorXd lambda0;
  get_lambda(V, o, lambda0);

  // L_tilda0
  Eigen::SparseMatrix<double> L_tilda0;
  laplacian(V, F, L_tilda0);

  // D_v
  Eigen::DiagonalMatrix<double, Eigen::Dynamic> D_v;
  sparse_v(V, D_v);

  // L_theta
  Eigen::VectorXd S_lambda0 = S * lambda0;
  Eigen::VectorXd L_theta;
  get_L_theta(S_lambda0, L_tilda0, D_v, L_theta);
  assert(L_theta.size() == 3 * V.rows());

  // D_L_theta
  Eigen::DiagonalMatrix<double, Eigen::Dynamic> D_L_theta_diag;
  D_L_theta_diag.diagonal() = L_theta;
  Eigen::SparseMatrix<double> D_L_theta;
  D_L_theta = D_L_theta_diag;

  Eigen::SparseMatrix<double> A;
  A = D_A * (L_tilda0 * D_v - D_L_theta) * S;

}

void smooth(
    const Eigen::MatrixXd & V,
    const Eigen::MatrixXi & F,
    const Eigen::MatrixXd & G,
    Eigen::MatrixXd & U)
{
  U.resize(G.rows(), G.cols());
}
