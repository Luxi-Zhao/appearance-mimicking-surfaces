#include "smooth.h"
#include <igl/read_triangle_mesh.h>
#include <igl/parula.h>
#include <igl/opengl/glfw/Viewer.h>
#include <igl/ambient_occlusion.h>
#include <igl/unproject_onto_mesh.h>
#include <string>
#include <iostream>
#include <igl/edge_lengths.h>
#include <igl/find.h>
#include <igl/cotmatrix.h>
#include <igl/repdiag.h>
#include "cotmatrix.h"
#include "massmatrix.h"

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
    (argc>1?argv[1]:"../data/cube.obj"),V,F);
  igl::opengl::glfw::Viewer viewer;
  std::cout<<R"(
[space]  Toggle whether displaying 3D surface or 2D parameterization
C,c      Toggle checkerboard
t        Switch parameterization to Tutte embedding
l        Switch parameterization to Least squares conformal mapping
)";

  const Eigen::RowVector3d orange(1.0,0.7,0.2);
  const Eigen::RowVector3d yellow(1.0,0.9,0.2);
  const Eigen::RowVector3d blue(0.2,0.3,0.8);
  const Eigen::RowVector3d green(0.2,0.6,0.3);

  // Find the bounding box
  Eigen::Vector3d m = V.colwise().minCoeff();
  Eigen::Vector3d M = V.colwise().maxCoeff();

  // Predetermined view points based on bounding box
  // TODO add the other 7 points
  Eigen::MatrixXd V_vp(1,3);
  V_vp <<
       (M(0)+m(0))/2.0, (M(1)+m(1))/2.0, 10.0 * M(2);

  // Lambda bounds according to view point
  // TODO add the other 7 points
  Eigen::MatrixXd lambda_bounds(1, 2);
  double lambda_hi = (V_vp.row(0) - Eigen::RowVector3d(m(0), m(1), m(2))).norm();
  double lambda_lo = lambda_hi * 0.93;

  bool deforming = false;

  const auto & update = [&]()
  {
      if(deforming)
      {
        deform(V, F, V_vp.row(0), lambda_lo, lambda_hi, DV);
        viewer.data().set_vertices(DV);
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
  Eigen::MatrixXd V_box(9,3);
  V_box <<
    m(0), m(1), m(2),
    M(0), m(1), m(2),
    M(0), M(1), m(2),
    m(0), M(1), m(2),
    m(0), m(1), M(2),
    M(0), m(1), M(2),
    M(0), M(1), M(2),
    m(0), M(1), M(2),
    V_vp.row(0);

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
  viewer.data().set_colors(yellow);
  update();
  viewer.data().show_texture = true;
  viewer.data().show_lines = false;
  viewer.launch();

  return EXIT_SUCCESS;
}
