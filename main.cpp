#include "smooth.h"
#include <igl/read_triangle_mesh.h>
#include <igl/write_triangle_mesh.h>
#include <igl/parula.h>
#include <igl/opengl/glfw/Viewer.h>
#include <igl/ambient_occlusion.h>
#include <string>
#include <iostream>

void get_weights(
  const Eigen::MatrixXd & V,
  const Eigen::MatrixXi & F,
  const Eigen::DiagonalMatrix<double, Eigen::Dynamic> & D_w)
{
  Eigen::MatrixXd N;
  igl::per_vertex_normals(V,F,N);
  Eigen::VectorXd S; // between 1 (fully occluded) and 0 (not occluded)
  igl::ambient_occlusion(V, F, V, N, 500, S);
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
[space]  Toggle whether displaying 3D surface or 2D parameterization
C,c      Toggle checkerboard
t        Switch parameterization to Tutte embedding
l        Switch parameterization to Least squares conformal mapping
)";

  const Eigen::RowVector3d white(1.0,1.0,1.0);

  // Find the bounding box
  Eigen::Vector3d m = V.colwise().minCoeff();
  Eigen::Vector3d M = V.colwise().maxCoeff();

  // Predetermined view point
  Eigen::RowVector3d o((M(0)+m(0))/2.0, (M(1)+m(1))/2.0, 2.0 * M(2));

  int dimen = 2; // z-axis
  double v_fixed = V(0, dimen); // fix the 0th vertex
  // vertex range along the z-axis
  double v_min = m(dimen);
  double v_max = M(dimen);

  double orig_depth = v_max - v_min;
  double half_depth = orig_depth / 6.0;
  double mesh_midpoint = (v_max + v_min) / 2.0;
  double lambda_max = mesh_midpoint - half_depth;
  double lambda_min = mesh_midpoint + half_depth;

  Eigen::MatrixXd V_hat = (V.rowwise() - o).rowwise().normalized();

  double pos_fixed = lambda_max + (v_fixed - v_min) / (v_max - v_min) * (lambda_min - lambda_max);
  double lambda_known = (pos_fixed - o(2)) / V_hat(0, 2);

  Eigen::VectorXd lambda_lo, lambda_hi;
  lambda_lo = ((lambda_min - o(2)) / V_hat.col(2).array()).matrix(); // TODO generalize this
  lambda_hi = ((lambda_max - o(2)) / V_hat.col(2).array()).matrix(); // TODO generalize this

  // TODO clean up here
  std::cout << "lambda known: " << lambda_known << "\n";
  std::cout << "lambda hi: " << lambda_hi[0] << "\n";
  std::cout << "lambda lo: " << lambda_lo[0] << "\n";
  if(lambda_known <= lambda_hi[0] && lambda_known >= lambda_lo[0]) {
    std::cout << "lambda known is good" << "\n";
  }

  bool deforming = false;

  const auto & update = [&]()
  {
      if(deforming)
      {
        deform(V, F, o, lambda_lo, lambda_hi, lambda_known, DV);
        viewer.data().set_vertices(DV);
        igl::write_triangle_mesh("../data/output.obj",DV,F);
      } else
      {
        viewer.data().set_vertices(V);
      }
      viewer.data().compute_normals();

  };
  viewer.callback_key_pressed =
    [&](igl::opengl::glfw::Viewer &, unsigned int key, int)
    {
        switch(key)
        {
          case ' ':
            deforming ^= 1;
            break;
          case 'C':
          case 'c':
            viewer.data().show_texture ^= 1;
            break;
          case 't':
            std::cout << viewer.core().camera_eye << "\n";
            std::cout << viewer.core().light_position << "\n";
            break;
          default:
            return false;
        }
        update();
        return true;
    };


  // Corners of the bounding box
  // TODO change this back to 8 later, vp is for testing only
  Eigen::MatrixXd V_box(10,3);
  V_box <<
    m(0), m(1), m(2),
    M(0), m(1), m(2),
    M(0), M(1), m(2),
    m(0), M(1), m(2),
    m(0), m(1), M(2),
    M(0), m(1), M(2),
    M(0), M(1), M(2),
    m(0), M(1), M(2),
    o,
    V.row(0);

  // Edges of the bounding box
  // TODO change this back to 12 later
  Eigen::MatrixXi E_box(13,2);
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
    7 ,3,
    0, 8;

  // Plot the corners of the bounding box as points
  viewer.data().add_points(V_box,Eigen::RowVector3d(1,0,0));

  // Plot the edges of the bounding box
  for (unsigned i=0;i<E_box.rows(); ++i)
    viewer.data().add_edges
      (
        V_box.row(E_box(i,0)),
        V_box.row(E_box(i,1)),
        Eigen::RowVector3d(1,0,0)
      );


  viewer.data().set_mesh(V,F);
  Eigen::MatrixXd N;
  igl::per_vertex_normals(V,F,N);
  viewer.data().set_colors(white);
  update();
  viewer.data().show_texture = true;
  viewer.data().show_lines = false;
  viewer.launch();

  return EXIT_SUCCESS;
}
