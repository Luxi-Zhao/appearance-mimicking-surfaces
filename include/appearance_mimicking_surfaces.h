#ifndef SMOOTH_H
#define SMOOTH_H
#include <Eigen/Core>
#include <Eigen/Sparse>

// Given a mesh, a viewpoint, and linear constraints:
//      lambda_lo -- combined with lambdaMax this will "squeeze" the mesh into the defined thickness by projecting the vertices
//      lambda_hi
//      mu_ind -- different parts of the shape can be given different thickness by segmenting the vertices into groups
//      fixed vertices
//      weights -- scales the difference of vertex normals so that more visible vertices are given preserved "more" than less visible
//
// output a bas-relief deformed mesh whose appearance is preserved
// from the viewpoint
//
//
// Inputs:
//    V            #V by 3 list of the vertex positions of the model
//    F            #F by 3 list of triangle indices into V
//    o            3D vector of the coordinates of the viewpoint
//    lambda_lo    #V by 1 list of lambdaMinValue, one for each vertex
//    lambda_hi    #V by 1 list of lambdaMaxValue, one for each vertex
//    ind_fixed    Index of the fixed vertex
//    lambda_known Lambda value of the fixed vertex
//    weights      #V length list of vertex weights
//    mu_ind       #V length list of mu indices
//                      (separates mesh into independent regions each with their own thickness constraint)
//                      (mu_ind[0] = 0 means the 0th vertex belongs to group 0)
//                      (default are all 0's, meaning there is only 1 group)
//  Outputs:
//    DV           #V by 3 list of the vertex positions of the deformed model

void appearance_mimicking_surfaces(
  const Eigen::MatrixXd & V,
  const Eigen::MatrixXi & F,
  const Eigen::RowVector3d & o,
  const Eigen::VectorXd & lambda_lo,
  const Eigen::VectorXd & lambda_hi,
  const Eigen::VectorXi & ind_fixed,
  double lambda_known,
  const Eigen::VectorXd & weights,
  const Eigen::VectorXi & mu_ind,
  Eigen::MatrixXd & DV);
#endif
