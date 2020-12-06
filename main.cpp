#include "appearance_mimicking_surfaces.h"
#include <igl/read_triangle_mesh.h>
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

void get_viewpoint_lower_mid(
  const Eigen::Vector3d & m,
  const Eigen::Vector3d & M,
  Eigen::RowVector3d & o)
{
  o = Eigen::RowVector3d((M(0)+m(0))/2.0, (M(1)+m(1))/2.0, 3.0 * M(2));
}

void get_viewpoint_lower_left(
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

void get_viewpoint_lower_right(
  const Eigen::Vector3d & m,
  const Eigen::Vector3d & M,
  Eigen::RowVector3d & o)
{
  Eigen::RowVector3d lower_right, upper_left, diagonal, o_low;
  lower_right = Eigen::RowVector3d(M(0), m(1), M(2));
  upper_left = m;
  diagonal = lower_right - upper_left;
  o_low = upper_left + diagonal * 2.0;
  o << o_low(0), (M(1)+m(1))/2.0, o_low(2);
}

void get_plane_mid(
  const Eigen::Vector3d & m,
  const Eigen::Vector3d & M,
  Eigen::Matrix3d & pl)
{
  Eigen::RowVector3d left_lo, right_lo, left_hi;
  left_lo = Eigen::RowVector3d(m(0), m(1), (m(2) + M(2)) / 2.0);
  right_lo = Eigen::RowVector3d(M(0), m(1), (m(2) + M(2)) / 2.0);
  left_hi = Eigen::RowVector3d(m(0), M(1), (m(2) + M(2)) / 2.0);
  pl.resize(3, 3);
  pl << left_hi, left_lo, right_lo;
}

void get_plane_left(
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

void get_plane_right(
  const Eigen::Vector3d & m,
  const Eigen::Vector3d & M,
  Eigen::Matrix3d & pl)
{
  Eigen::RowVector3d lower_left_lo, upper_right_lo, lower_left_hi;
  lower_left_lo = Eigen::RowVector3d(m(0), m(1), M(2));
  upper_right_lo = Eigen::RowVector3d(M(0), m(1), m(2));
  lower_left_hi = Eigen::RowVector3d(m(0), M(1), M(2));
  pl.resize(3, 3);
  pl << lower_left_hi, lower_left_lo, upper_right_lo;
}

double get_bounds(
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

  double lambda_known = (lambda_lo(ind_fixed) + lambda_hi(ind_fixed)) / 2.0;
  return lambda_known;
}

int main(int argc, char *argv[])
{
  // Load input meshes
  Eigen::MatrixXd V, DV_m, DV_l, DV_r;
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

  int ind_fixed = 0;
  Eigen::RowVector3d o, o_m, o_l, o_r;
  Eigen::Matrix3d pl_m, pl_l, pl_r;
  Eigen::VectorXd lambda_lo_m, lambda_lo_l, lambda_lo_r, lambda_hi_m, lambda_hi_l, lambda_hi_r;
  double lambda_known_m, lambda_known_l, lambda_known_r;

  get_viewpoint_lower_mid(m, M, o_m);
  get_viewpoint_lower_left(m, M, o_l);
  get_viewpoint_lower_right(m, M, o_r);
  o = o_m;

  get_plane_mid(m, M, pl_m);
  get_plane_left(m, M, pl_l);
  get_plane_right(m, M, pl_r);

  lambda_known_m = get_bounds(o_m, V, pl_m, ind_fixed, lambda_lo_m, lambda_hi_m);
  lambda_known_l = get_bounds(o_l, V, pl_l, ind_fixed, lambda_lo_l, lambda_hi_l);
  lambda_known_r = get_bounds(o_r, V, pl_r, ind_fixed, lambda_lo_r, lambda_hi_r);

  Eigen::VectorXd weights(V.rows());
  weights.setOnes();

  Eigen::VectorXi mu_ind(V.rows());
  mu_ind.setZero();

  appearance_mimicking_surfaces(V, F, o_m, lambda_lo_m, lambda_hi_m, ind_fixed, lambda_known_m, weights, mu_ind, DV_m);
  appearance_mimicking_surfaces(V, F, o_l, lambda_lo_l, lambda_hi_l, ind_fixed, lambda_known_l, weights, mu_ind, DV_l);
  appearance_mimicking_surfaces(V, F, o_r, lambda_lo_r, lambda_hi_r, ind_fixed, lambda_known_r, weights, mu_ind, DV_r);

  // 0 - do not show, 1 - show mid view,
  // 2 - show left view, 3 - show right view
  int show_deform = 0;
  bool show_points = false;

  const auto & update = [&]()
  {
      if(show_points)
      {
        // Plot the corners of the bounding box as points
        viewer.data().add_points(V_corners, red);
        viewer.data().add_points(o, green);
        viewer.data().add_points(V.row(0), blue);
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

      if(show_deform == 1)
      {
        viewer.data().set_vertices(DV_m);
      } else if(show_deform == 2)
      {
        viewer.data().set_vertices(DV_l);
      } else if(show_deform == 3)
      {
        viewer.data().set_vertices(DV_r);
      } else {
        viewer.data().set_vertices(V);
      }
  };
  viewer.callback_key_pressed =
    [&](igl::opengl::glfw::Viewer &, unsigned int key, int)
    {
        switch(key)
        {
          case 'm':
            show_deform = 1;
            o = o_m;
            break;
          case 'l':
            show_deform = 2;
            o = o_l;
            break;
          case 'r':
            show_deform = 3;
            o = o_r;
            break;
          case ' ':
            show_deform = 0;
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
