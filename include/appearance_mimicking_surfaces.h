#ifndef SMOOTH_H
#define SMOOTH_H
#include <Eigen/Core>
#include <Eigen/Sparse>

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
  Eigen::MatrixXd & DV);
#endif
