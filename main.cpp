#include "appearance_mimicking_surfaces.h"
#include <igl/read_triangle_mesh.h>
#include <igl/write_triangle_mesh.h>
#include <igl/parula.h>
#include <igl/opengl/glfw/Viewer.h>
#include <string>
#include <iostream>

double point_plane_distance(
  const Eigen::RowVector3d & x,
  const Eigen::RowVector3d & a,
  const Eigen::RowVector3d & b,
  const Eigen::RowVector3d & c,
  Eigen::RowVector3d & p)
{
  // Project point x onto the plane pl containing triangle abc
  Eigen::RowVector3d n_p;
  n_p = (b-a).cross(c-a);
  double d_xpl = (x-a).dot(n_p) / n_p.norm();
  p = x - n_p.normalized() * d_xpl;
  double d = (x-p).norm();
  return d;
}

double mesh_depth(
  const Eigen::MatrixXd & V,
  int dimen)
{
  // Vertex range along the z-axis
  double v_min = V.col(dimen).minCoeff();
  double v_max = V.col(dimen).maxCoeff();
  double d = v_max - v_min;
  return d;
}

void get_viewpoint(
  const Eigen::Vector3d & m,
  const Eigen::Vector3d & M,
  Eigen::RowVector3d & o)
{
  Eigen::RowVector3d lower_left, upper_right, diagonal, o_low;
  lower_left = Eigen::RowVector3d(m(0), m(1), M(2));
  upper_right = Eigen::RowVector3d(M(0), m(1), m(2));
  diagonal = lower_left - upper_right;
  o_low = upper_right + diagonal * 2.0;
  o << o_low(0), (M(1)+m(1))/2.0, o_low(2);
}

void get_plane(
  const Eigen::Vector3d & m,
  const Eigen::Vector3d & M,
  Eigen::Matrix3d & pl)
{
  Eigen::RowVector3d upper_left_lo, lower_right_lo, lower_right_hi;
  upper_left_lo = m;
  lower_right_lo = Eigen::RowVector3d(M(0), m(1), M(2));
  lower_right_hi = Eigen::RowVector3d(M(0), M(1), M(2));
  pl.resize(3, 3);
  pl << lower_right_hi, upper_left_lo, lower_right_lo;
}

double get_bounds2(
  const Eigen::RowVector3d & o,
  const Eigen::MatrixXd & V,
  const Eigen::Matrix3d & pl,
  int ind_fixed,
  Eigen::VectorXd & lambda_lo,
  Eigen::VectorXd & lambda_hi)
{
  Eigen::MatrixXd V_hat = (V.rowwise() - o).rowwise().normalized();
  double orig_depth = mesh_depth(V, 2);
  double scale = 1.0/3.0;
  double half_depth = scale * orig_depth / 2.0;

  Eigen::RowVector3d a, b, c, p;
  a = pl.row(0);
  b = pl.row(1);
  c = pl.row(2);
  double mid_dist = point_plane_distance(o, a, b, c, p);
  double lambda_max = mid_dist + half_depth;
  double lambda_min = mid_dist - half_depth;

  Eigen::RowVector3d op_hat = (p - o).normalized();
  Eigen::VectorXd denom(V_hat.rows());
  for(int i = 0; i < V_hat.rows(); i++) {
    denom(i) = V_hat.row(i).dot(op_hat);
  }

  lambda_lo = lambda_min / denom.array();
  lambda_hi = lambda_max / denom.array();

  // TODO properly calculate this
  double lambda_known = (lambda_min + lambda_max) / 2.0;
  return lambda_known;
}


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

  // Colors
  const Eigen::RowVector3d black(0.0,0.0,0.0);
  const Eigen::RowVector3d white(1.0,1.0,1.0);
  const Eigen::RowVector3d red(1.0,0.0,0.0);
  const Eigen::RowVector3d green(0.0,1.0,0.0);
  const Eigen::RowVector3d blue(0.0,0.0,1.0);

  // Find a bounding box for the mesh
  Eigen::Vector3d m = V.colwise().minCoeff();
  Eigen::Vector3d M = V.colwise().maxCoeff();

  Eigen::MatrixXd V_corners(2,3);
  V_corners << m(0), m(1), m(2),
    M(0), M(1), M(2);

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

  // Position view point along the z-axis
//  Eigen::RowVector3d o((M(0)+m(0))/2.0, (M(1)+m(1))/2.0, 2.0 * M(2));
  Eigen::RowVector3d o;
  get_viewpoint(m, M, o);
  Eigen::Matrix3d pl;
  get_plane(m, M, pl);

  int ind_fixed = 0;
  Eigen::VectorXd lambda_lo, lambda_hi, weights(V.rows());
  double lambda_known = get_bounds2(o, V, pl, ind_fixed, lambda_lo, lambda_hi);
  weights.setOnes();

  Eigen::VectorXi mu_ind(V.rows());
  mu_ind.setZero();

  appearance_mimicking_surfaces(V, F, o, lambda_lo, lambda_hi, ind_fixed, lambda_known, weights, mu_ind, DV);
//  igl::write_triangle_mesh("../data/output.obj",DV,F);

  bool show_deform = false;
  bool show_points = false;

  const auto & update = [&]()
  {
      if(show_points)
      {
        // Plot the corners of the bounding box as points
//        viewer.data().add_points(V_box, red);
        viewer.data().add_points(o, green);
        viewer.data().add_points(V.row(0), blue);
        viewer.data().add_points(V_corners, black);
        // Plot the edges of the bounding box
        for (unsigned i = 0; i < E_box.rows(); ++i)
          viewer.data().add_edges
            (
              V_box.row(E_box(i,0)),
              V_box.row(E_box(i,1)),
              red
            );
      } else {
        viewer.data().clear_points();
        viewer.data().clear_edges();
      }

      if(show_deform)
      {
        viewer.data().set_vertices(DV);
      } else
      {
        viewer.data().set_vertices(V);
      }
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
          default:
            return false;
        }
        update();
        return true;
    };

  viewer.data().set_mesh(V,F);
  viewer.data().set_colors(white);
  viewer.core().background_color.setOnes();
  update();
  viewer.data().show_texture = true;
  viewer.data().show_lines = false;
  viewer.launch();

  return EXIT_SUCCESS;
}
