#include "smooth.h"
#include <igl/cotmatrix.h>
#include <igl/massmatrix.h>
#include <igl/repdiag.h>
#include <igl/invert_diag.h>
#include <igl/min_quad_with_fixed.h>
#include <igl/active_set.h>
#include <cmath>
#include <iostream>

typedef Eigen::Triplet<double> T;

// D_A
void voronoi_area(
  const Eigen::SparseMatrix<double> & M,
  Eigen::DiagonalMatrix<double, Eigen::Dynamic> & D_A)
{
  Eigen::VectorXd M_diag = M.diagonal();

  D_A.resize(M_diag.size() * 3);
  D_A.setZero();

  for(int i = 0; i < M_diag.size(); i++) {
    for(int j = 0; j < 3; j++) {
      D_A.diagonal()[3 * i + j] = sqrt(M_diag[i]);
    }
  }
}

void get_weights(
  const Eigen::VectorXd & weights,
  Eigen::DiagonalMatrix<double, Eigen::Dynamic> & D_w)
{
  D_w.resize(weights.size() * 3);
  D_w.setZero();
  for(int i = 0; i < weights.size(); i++) {
    for(int j = 0; j < 3; j++) {
      D_w.diagonal()[3 * i + j] = weights(i);
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

void diag_concat(
  const Eigen::SparseMatrix<double> & A,
  const Eigen::SparseMatrix<double> & B,
  Eigen::SparseMatrix<double> & C)
{
  Eigen::SparseMatrix<double> zeros_A, zeros_B, row1, row2;
  zeros_A.resize(A.rows(), B.cols());
  zeros_B.resize(B.rows(), A.cols());

  igl::cat(2, A, zeros_A, row1);
  igl::cat(2, zeros_B, B, row2);
  igl::cat(1, row1, row2, C);
}


void deform(
  const Eigen::MatrixXd & V,
  const Eigen::MatrixXi & F,
  const Eigen::RowVector3d & o,
  const Eigen::VectorXd & lambda_lo,
  const Eigen::VectorXd & lambda_hi,
  int ind_fixed,
  double lambda_known,
  const Eigen::VectorXd & weights,
  const Eigen::VectorXi & mu_ind,
  Eigen::MatrixXd & DV)
{
  double test = (lambda_hi.array()-lambda_lo.array()).minCoeff();
  if(test <= 0) {
    std::cout << "ASSERT FAILED!!!!!!!!" << std::endl;
  } else {
    std::cout << "assert passed" << std::endl;
  }
  assert((lambda_hi.array()-lambda_lo.array()).minCoeff() > 0 && "ux(i) must be > lx(i)");

  Eigen::SparseMatrix<double> M;
  igl::massmatrix(V, F, igl::MASSMATRIX_TYPE_VORONOI, M);

  std::cout << "getting D_A" << std::endl;
  // D_A
  Eigen::DiagonalMatrix<double, Eigen::Dynamic> D_A;
  voronoi_area(M, D_A);

  std::cout << "getting D_w" << std::endl;

  // D_w
  Eigen::DiagonalMatrix<double, Eigen::Dynamic> D_w;
  get_weights(weights, D_w);
  Eigen::SparseMatrix<double> D_W = D_w.toDenseMatrix().sparseView();

  std::cout << "getting S" << std::endl;
  // S
  Eigen::SparseMatrix<double> S;
  selector(V, S);

  std::cout << "getting lambda0" << std::endl;

  // lambda0
  Eigen::VectorXd lambda0;
  get_lambda(V, o, lambda0);

  std::cout << "getting L_tilda0" << std::endl;

  // L_tilda0
  Eigen::SparseMatrix<double> L_tilda0;
  laplacian(V, F, M, L_tilda0);

  std::cout << "getting D_v_hat" << std::endl;

  // D_v_hat
  Eigen::DiagonalMatrix<double, Eigen::Dynamic> D_v_hat;
  Eigen::MatrixXd V_hat = (V.rowwise() - o).rowwise().normalized();
  sparse_v(V_hat, D_v_hat);

  std::cout << "getting L_theta" << std::endl;

  // L_theta
  Eigen::VectorXd S_lambda0 = S * lambda0;
  Eigen::VectorXd L_theta;
  get_L_theta(S_lambda0, L_tilda0, D_v_hat, L_theta);

  std::cout << "getting D_L_theta" << std::endl;

  // D_L_theta
  Eigen::DiagonalMatrix<double, Eigen::Dynamic> D_L_theta_diag;
  D_L_theta_diag.diagonal() = L_theta;
  Eigen::SparseMatrix<double> D_L_theta;
  D_L_theta = D_L_theta_diag;

  std::cout << "getting Q" << std::endl;

  Eigen::SparseMatrix<double> Q;
  Q = D_A * D_W * (L_tilda0 * D_v_hat - D_L_theta) * S;

  std::cout << "getting F" << std::endl;

  // F
  int n_mu = mu_ind.maxCoeff() + 1;
  double alpha_sqrt = sqrt(pow(10.0, -7.0));
  Eigen::SparseMatrix<double> I = alpha_sqrt * Eigen::MatrixXd::Identity(n_mu, n_mu).sparseView();
  Eigen::SparseMatrix<double> J, K;
  diag_concat(Q, I, J);
  K = J.transpose() * J;

  std::cout << "preparing constraints" << std::endl;
  long n = V.rows() + n_mu;
  assert(K.rows() == n);
  assert(K.cols() == n);
  // There are no linear coefficients so set B to 0
  Eigen::MatrixXd B = Eigen::MatrixXd::Zero(n, 1);

  // Fix the value for lambda for one vertex
  // to get a unique solution
  Eigen::VectorXi b(1);
  Eigen::VectorXd Y(1);
  b(0) = ind_fixed;
  Y(0) = lambda_known;

  Eigen::SparseMatrix<double> Aeq, Aieq;
  Eigen::VectorXd Beq, Bieq;
  Eigen::VectorXd lx(n), ux(n), mu_l, mu_u;

  // Not really sure what the bounds should
  // be for mu's
  mu_l = Eigen::VectorXd::Ones(n_mu) * 1.0;
  mu_u = Eigen::VectorXd::Ones(n_mu) * 1.5;
  lx << lambda_lo, mu_l;
  ux << lambda_hi, mu_u;

  igl::active_set_params as;
  Eigen::VectorXd x;

  std::cout << "solving for x's" << std::endl;
  igl::active_set(K, B, b, Y, Aeq, Beq, Aieq, Bieq, lx, ux, as, x);

  std::cout << "got x's" << std::endl;
  std::cout << x << std::endl;

  // Extract lambda and mu from x
  Eigen::VectorXd lambda, mu, mu_arr(V.rows());
  lambda = x.head(V.rows());
  mu = x.tail(n_mu);
  for(int i = 0; i < mu_arr.size(); i++) {
    mu_arr(i) = mu(mu_ind(i));
  }

  // Use lambdas to transform vertices
  std::cout << "getting DV" << std::endl;
  DV.resize(V.rows(), V.cols());
  DV = (V_hat.array().colwise() * (lambda.array() * mu_arr.array())).matrix();
  DV = DV.rowwise() + o;

  std::cout << "DONE" << std::endl;
}

void smooth(
    const Eigen::MatrixXd & V,
    const Eigen::MatrixXi & F,
    const Eigen::MatrixXd & G,
    Eigen::MatrixXd & U)
{
  U.resize(G.rows(), G.cols());
}
