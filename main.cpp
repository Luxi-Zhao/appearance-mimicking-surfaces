#include "appearance_mimicking_surfaces.h"
#include <igl/read_triangle_mesh.h>
#include <igl/write_triangle_mesh.h>
#include <igl/parula.h>
#include <igl/png/readPNG.h>
#include <igl/opengl/glfw/Viewer.h>
#include <string>
#include <iostream>

/**
 * Compute lambda constraints according to view point position.
 * This method is not generalized for all meshes. May need to
 * adjust it for different meshes.
 */
double get_bounds(
  const Eigen::RowVector3d & o,
  const Eigen::MatrixXd & V,
  int ind_fixed,
  Eigen::VectorXd & lambda_lo,
  Eigen::VectorXd & lambda_hi)
{
  int dimen = 2; // z-axis

  Eigen::MatrixXd V_hat = (V.rowwise() - o).rowwise().normalized();
  double v_fixed = V(ind_fixed, dimen);

  // Vertex range along the z-axis
  double v_min = V.col(dimen).minCoeff();
  double v_max = V.col(dimen).maxCoeff();
  double orig_depth = v_max - v_min;

  double scale = 1.0/3.0;
  double half_depth = scale * orig_depth / 2.0;
  double mesh_midpoint = (v_max + v_min) / 2.0;
  // Constrain the deformed mesh between planes defined by
  // these two bounds
  double lambda_max = mesh_midpoint - half_depth;
  double lambda_min = mesh_midpoint + half_depth;

  // New position of the 0th vertex after deformation
  double pos_fixed = lambda_max + (v_fixed - v_min) / (v_max - v_min) * (lambda_min - lambda_max);
  double lambda_known = (pos_fixed - o(dimen)) / V_hat(ind_fixed, dimen);

  lambda_lo = ((lambda_min - o(dimen)) / V_hat.col(dimen).array()).matrix();
  lambda_hi = ((lambda_max - o(dimen)) / V_hat.col(dimen).array()).matrix();

  return lambda_known;
}

void get_background(
  const Eigen::MatrixXd & V,
  const Eigen::MatrixXi & F,
  Eigen::MatrixXd & MV,
  Eigen::MatrixXi & MF)
{
  Eigen::Vector3d m = V.colwise().minCoeff();
  Eigen::Vector3d M = V.colwise().maxCoeff();

  // Corners of background
  Eigen::MatrixXd V_back(4,3);
  double half_x = (M(0) - m(0)) / 2.0;
  double half_y = (M(1) - m(1)) / 2.0;
  V_back <<
         m(0) - half_x, m(1) - half_y, m(2),
    M(0) + half_x, m(1) - half_y, m(2),
    M(0) + half_x, M(1) + half_y, m(2),
    m(0) - half_x, M(1) + half_y, m(2);
  int ind_start = V.rows();
  Eigen::MatrixXi F_back(2, 3);
  F_back <<
         ind_start, ind_start + 1, ind_start + 2,
    ind_start, ind_start + 2, ind_start + 3;

  MV.resize(V.rows() + 4, V.cols());
  MF.resize(F.rows() + 2, F.cols());
  MV << V, V_back;
  MF << F, F_back;
}

int main(int argc, char *argv[])
{
  // Load input meshes
  Eigen::MatrixXd V,DV;
  Eigen::MatrixXi F;
  igl::read_triangle_mesh(
    (argc>1?argv[1]:"../data/icosphere.obj"),V,F);
  igl::opengl::glfw::Viewer viewer;
  std::cout<<R"(
[space]  Toggle whether displaying original mesh or deformed mesh
p        Toggle debug points (red - bounding box, green - view point, blue - fixed vertex)
)";
  Eigen::Matrix<unsigned char,Eigen::Dynamic,Eigen::Dynamic> R,G,B,A;
  igl::png::readPNG(argv[2],R,G,B,A);

  // Colors
  const Eigen::RowVector3d white(1.0,1.0,1.0);
  const Eigen::RowVector3d red(1.0,0.0,0.0);
  const Eigen::RowVector3d green(0.0,1.0,0.0);
  const Eigen::RowVector3d blue(0.0,0.0,1.0);

  // Find a bounding box for the mesh
  Eigen::Vector3d m = V.colwise().minCoeff();
  Eigen::Vector3d M = V.colwise().maxCoeff();

  // Position view point along the z-axis
  Eigen::RowVector3d o((M(0)+m(0))/2.0, (M(1)+m(1))/2.0, 2.0 * M(2));

  // Corners of the bounding box
  Eigen::MatrixXd V_box(8,3);
  V_box <<
    m(0), m(1), m(2),
    M(0), m(1), m(2),
    M(0), M(1), m(2),
    m(0), M(1), m(2),
    m(0), m(1), M(2),
    M(0), m(1), M(2),
    M(0), M(1), M(2),
    m(0), M(1), M(2);

  // Edges of the bounding box
  Eigen::MatrixXi E_box(12,2);
  E_box <<
    0, 1,
    1, 2,
    2, 3,
    3, 0,
    4, 5,
    5, 6,
    6, 7,
    7, 4,
    0, 4,
    1, 5,
    2, 6,
    7 ,3;

  int ind_fixed = 0;
  Eigen::VectorXd lambda_lo, lambda_hi, weights(V.rows());
  double lambda_known = get_bounds(o, V, ind_fixed, lambda_lo, lambda_hi);
  weights.setOnes();

  Eigen::VectorXi mu_ind(V.rows());
  mu_ind.setZero();

//  deform(V, F, o, lambda_lo, lambda_hi, ind_fixed, lambda_known, weights, mu_ind, DV);
//  igl::write_triangle_mesh("../data/output.obj",DV,F);

  bool show_deform = false;
  bool show_points = false;
  bool show_back = false;

  Eigen::MatrixXd CV;
  CV = V;

  const auto & update = [&]()
  {
      viewer.data().clear();
      Eigen::MatrixXd MV;
      Eigen::MatrixXi MF;
      if(show_deform && show_back) {
        get_background(DV, F, MV, MF);
        viewer.data().set_mesh(MV, MF);
      } else if(show_deform){
        viewer.data().set_mesh(DV, F);
      } else if(show_back) {
        get_background(V, F, MV, MF);
        viewer.data().set_mesh(MV, MF);
      } else {
        viewer.data().set_mesh(V, F);
      }

      if(show_points) {
        // Plot the corners of the bounding box as points
        viewer.data().add_points(V_box, red);
        viewer.data().add_points(o, green);
        viewer.data().add_points(V.row(0), blue);
        // Plot the edges of the bounding box
        for (unsigned i = 0; i < E_box.rows(); ++i)
          viewer.data().add_edges
            (
              V_box.row(E_box(i, 0)),
              V_box.row(E_box(i, 1)),
              red
            );
      }

      viewer.data().set_texture(R,G,B,A);
      viewer.data().use_matcap = true;
      viewer.data().face_based = true;
  };
  viewer.callback_key_pressed =
    [&](igl::opengl::glfw::Viewer &, unsigned int key, int)
    {
        switch(key)
        {
          case ' ':
            show_deform ^= 1;
            break;
          case 'p':
            show_points ^= 1;
            break;
          case 'b':
            show_back ^= 1;
            break;
          default:
            return false;
        }
        update();
        return true;
    };

  viewer.data().set_mesh(V,F);
  viewer.data().set_texture(R,G,B,A);
  viewer.data().use_matcap = true;
  viewer.core().background_color.setOnes();
  update();
  viewer.data().show_texture = true;
  viewer.data().show_lines = false;
  viewer.data().face_based = true;
  viewer.launch();

  return EXIT_SUCCESS;
}
