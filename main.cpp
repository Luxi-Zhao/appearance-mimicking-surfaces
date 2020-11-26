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

  //deform(V, F, DV);
  bool deform = false;
  bool place_viewpoint = false;

  Eigen::MatrixXd P(1,2);
//  P.row(0) = Eigen::RowVector3f(0, 0, 0);
//  P.row(1) = Eigen::RowVector3d(0, 0, 5);
  Eigen::RowVector3f last_mouse;

  const auto & update = [&]()
  {
      if(place_viewpoint) {
        viewer.data().set_points(P, yellow);
      }
      if(deform)
      {
        viewer.data().set_vertices(DV);
      }else
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
            deform ^= 1;
            break;
          case 'C':
          case 'c':
            viewer.data().show_texture ^= 1;
            break;
          case 't':
            std::cout << viewer.core().camera_eye << "\n";
            std::cout << viewer.core().light_position << "\n";
            break;
          case 'v':
            place_viewpoint ^= 1;
            break;
          default:
            return false;
        }
        update();
        return true;
    };

  viewer.callback_mouse_down =
    [&](igl::opengl::glfw::Viewer&, int, int)->bool
    {
        last_mouse = Eigen::RowVector3f(
          viewer.current_mouse_x,viewer.core().viewport(3)-viewer.current_mouse_y,0);
        if(place_viewpoint){
          // Find closest point on mesh to mouse position
          int fid;
          Eigen::Vector3f bary;
          if(igl::unproject_onto_mesh(
            last_mouse.head(2),
            viewer.core().view,
            viewer.core().proj,
            viewer.core().viewport,
            V, F,
            fid, bary))
          {
            long c;
            bary.maxCoeff(&c);
            Eigen::RowVector3d new_c = V.row(F(fid,c));
            if(P.size()==0 || (P.rowwise()-new_c).rowwise().norm().minCoeff() > 0)
            {
              P.conservativeResize(P.rows()+1,3);
              // Snap to closest vertex on hit face
              P.row(P.rows()-1) = new_c;
              update();
              return true;
            }
          }
        }
        return false;
    };

  viewer.data().set_mesh(V,F);
  Eigen::MatrixXd N;
  igl::per_vertex_normals(V,F,N);
  viewer.data().set_colors(N.array()*0.5+0.5);
  update();
  viewer.data().show_texture = true;
  viewer.data().show_lines = false;
  viewer.launch();

  return EXIT_SUCCESS;
}
