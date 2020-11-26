#include "smooth.h"
#include <igl/cotmatrix.h>
#include <igl/massmatrix.h>
#include <igl/repdiag.h>
#include <igl/invert_diag.h>
#include <igl/min_quad_with_fixed.h>
#include <igl/active_set.h>
#include <iostream>

typedef Eigen::Triplet<double> T;

// D_A
void voronoi_area(
  const Eigen::MatrixXd & V,
  const Eigen::MatrixXi & F,
  const Eigen::SparseMatrix<double> & M,
  Eigen::DiagonalMatrix<double, Eigen::Dynamic> & D_A)
{
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
  const Eigen::SparseMatrix<double> & M,
  Eigen::SparseMatrix<double> & L_tilda)
{
  Eigen::SparseMatrix<double> cot, M_inv, L;
  igl::cotmatrix(V, F, cot);
  igl::invert_diag(M, M_inv);
  L = M_inv * cot;

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

// D_v_hat
void sparse_v(
  const Eigen::MatrixXd & V,
  Eigen::DiagonalMatrix<double, Eigen::Dynamic> & D_v_hat)
{
  D_v_hat.resize(V.rows() * 3);
  D_v_hat.setZero();

  for(int i = 0; i < V.rows(); i++) {
    for(int j = 0; j < 3; j++) {
      D_v_hat.diagonal()[3 * i + j] = V(i, j);
    }
  }
}

// S
void selector(
  const Eigen::MatrixXd & V,
  Eigen::SparseMatrix<double> & S)
{
  Eigen::Vector3d ones(3);
  ones.setOnes();
  Eigen::SparseMatrix<double> ones_sp = ones.sparseView();

  S.resize(3 * V.rows(), V.rows());
  igl::repdiag(ones_sp, V.rows(), S);
}

void get_lambda(
  const Eigen::MatrixXd & V,
  const Eigen::RowVector3d & o,
  Eigen::VectorXd & lambda)
{
  lambda.resize(V.rows());
  lambda = (V.rowwise() - o).rowwise().norm();
}

void get_L_theta(
  const Eigen::VectorXd & S_lambda0,
  const Eigen::SparseMatrix<double> & L_tilda0,
  const Eigen::DiagonalMatrix<double, Eigen::Dynamic> & D_v_hat,
  Eigen::VectorXd & L_theta)
{
  Eigen::DiagonalMatrix<double, Eigen::Dynamic> D_S_lambda0;
  D_S_lambda0.diagonal() = S_lambda0;

  L_theta.setZero();
  L_theta = D_S_lambda0.inverse() * L_tilda0 * D_v_hat * S_lambda0;
}

void deform(
  const Eigen::MatrixXd & V,
  const Eigen::MatrixXi & F,
  const Eigen::RowVector3d & o,
  Eigen::MatrixXd & DV)
{
  std::cout << "testing deform" << std::endl;

  Eigen::SparseMatrix<double> M;
  igl::massmatrix(V, F, igl::MASSMATRIX_TYPE_VORONOI, M);

  std::cout << "getting D_A" << std::endl;
  // D_A
  Eigen::DiagonalMatrix<double, Eigen::Dynamic> D_A;
  voronoi_area(V, F, M, D_A);

  std::cout << "getting S" << std::endl;
  // S
  Eigen::SparseMatrix<double> S;
  selector(V, S);

  std::cout << "getting lambda0" << std::endl;

  // lambda0
//  Eigen::RowVector3d o;
//  o.setZero();
  Eigen::VectorXd lambda0;
  get_lambda(V, o, lambda0);
  std::cout << "---got lambda0----" << std::endl;
  std::cout << lambda0 << std::endl;

  std::cout << "getting L_tilda0" << std::endl;

  // L_tilda0
  Eigen::SparseMatrix<double> L_tilda0;
  laplacian(V, F, M, L_tilda0);

  std::cout << "getting D_v_hat" << std::endl;

  // D_v_hat
  Eigen::DiagonalMatrix<double, Eigen::Dynamic> D_v_hat;
  Eigen::MatrixXd V_hat = V.rowwise().normalized();
  sparse_v(V_hat, D_v_hat);

  std::cout << "getting L_theta" << std::endl;

  // L_theta
  Eigen::VectorXd S_lambda0 = S * lambda0;
  Eigen::VectorXd L_theta;
  get_L_theta(S_lambda0, L_tilda0, D_v_hat, L_theta);
  assert(L_theta.size() == 3 * V.rows());

  std::cout << "getting D_L_theta" << std::endl;

  // D_L_theta
  Eigen::DiagonalMatrix<double, Eigen::Dynamic> D_L_theta_diag;
  D_L_theta_diag.diagonal() = L_theta;
  Eigen::SparseMatrix<double> D_L_theta;
  D_L_theta = D_L_theta_diag;

  std::cout << "getting A" << std::endl;

  Eigen::SparseMatrix<double> A, Q;
  A = D_A * (L_tilda0 * D_v_hat - D_L_theta) * S;
  Q = A.transpose() * A;

  std::cout << "preparing constraints" << std::endl;
  // There are no linear coefficients so set B to 0
  Eigen::MatrixXd B = Eigen::MatrixXd::Zero(V.rows(), 1);

  double lambda_lo = 0.1;
  double lambda_hi = 0.2;
  // Fix the value for lambda for the 0th vertex
  Eigen::VectorXi b(1);
  Eigen::VectorXd Y(1);
  b(0) = 0;
  Y(0) = (lambda_lo + lambda_hi) / 2.0;

  Eigen::SparseMatrix<double> Aeq, Aieq;
  Eigen::VectorXd Beq, Bieq;
  Eigen::VectorXd lx, ux;
  lx = Eigen::VectorXd::Ones(V.rows(), 1) * lambda_lo;
  ux = Eigen::VectorXd::Ones(V.rows(), 1) * lambda_hi;

  igl::active_set_params as;
  Eigen::VectorXd lambda;

  std::cout << "solving for lambdas" << std::endl;
  igl::active_set(Q, B, b, Y, Aeq, Beq, Aieq, Bieq, lx, ux, as, lambda);

  std::cout << "got lambda" << std::endl;
  std::cout << lambda << std::endl;

  // Use lambdas to transform vertices
  std::cout << "getting DV" << std::endl;
  DV.resize(V.rows(), V.cols());
  DV = (V_hat.array().colwise() * lambda.array()).matrix();
}

void smooth(
    const Eigen::MatrixXd & V,
    const Eigen::MatrixXi & F,
    const Eigen::MatrixXd & G,
    Eigen::MatrixXd & U)
{
  U.resize(G.rows(), G.cols());
}
