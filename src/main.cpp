#include <igl/opengl/glfw/Viewer.h>

#include <igl/opengl/glfw/imgui/ImGuiPlugin.h>
#include <igl/opengl/glfw/imgui/ImGuiMenu.h>
#include <igl/opengl/glfw/imgui/ImGuiHelpers.h>

#include <igl/copyleft/tetgen/tetrahedralize.h>
#include <igl/readOBJ.h>
#include <igl/readPLY.h>
#include <igl/readSTL.h>
#include <igl/barycenter.h>
#include <iostream>

using namespace std;
using namespace Eigen;

#define DEBUG 0

// Input polygon
Eigen::MatrixXd V; // #V by 3 matrix for vertices
Eigen::MatrixXi F; // matrix for face indices
Eigen::MatrixXd B; // matrix for barycenters
Eigen::MatrixXd N; // matrix for normals

// Tetrahedralized interior
Eigen::MatrixXd TV; // #TV by 3 matrix for vertex positions
Eigen::MatrixXi TT; // #TT by 4 matrix for tet face indices
Eigen::MatrixXi TF; // #TF by 3 matrix for triangle face indices ('f', else `boundary_facets` is called on TT)

// Helper func, called every time a keyboard button is pressed
// Slice model at various percentage to view internal tetrahedral structure
bool key_down(igl::opengl::glfw::Viewer& viewer, unsigned char key, int modifier)
{
    if (key >= '1' && key <= '9')
    {
        double t = double((key - '1') + 1) / 9.0;

        VectorXd v = B.col(2).array() - B.col(2).minCoeff();
        v /= v.col(0).maxCoeff();

        vector<int> s;

        for (unsigned i = 0; i < v.size(); ++i)
            if (v(i) < t)
                s.push_back(i);

        MatrixXd V_temp(s.size() * 4, 3);
        MatrixXi F_temp(s.size() * 4, 3);
        MatrixXd B_temp(s.size(), 3);

        for (unsigned i = 0; i < s.size(); ++i)
        {
            V_temp.row(i * 4 + 0) = TV.row(TT(s[i], 0));
            V_temp.row(i * 4 + 1) = TV.row(TT(s[i], 1));
            V_temp.row(i * 4 + 2) = TV.row(TT(s[i], 2));
            V_temp.row(i * 4 + 3) = TV.row(TT(s[i], 3));
            F_temp.row(i * 4 + 0) << (i * 4) + 0, (i * 4) + 1, (i * 4) + 3;
            F_temp.row(i * 4 + 1) << (i * 4) + 0, (i * 4) + 2, (i * 4) + 1;
            F_temp.row(i * 4 + 2) << (i * 4) + 3, (i * 4) + 2, (i * 4) + 0;
            F_temp.row(i * 4 + 3) << (i * 4) + 1, (i * 4) + 2, (i * 4) + 3;
            B_temp.row(i) = (TV.row(TT(s[i], 0)) +
                TV.row(TT(s[i], 1)) +
                TV.row(TT(s[i], 2)) +
                TV.row(TT(s[i], 3))) / 4;
        }

        viewer.data().clear();
        viewer.data().set_mesh(V_temp, F_temp);
        viewer.data().add_points(B_temp, Eigen::RowVector3d(1, 0, 0));
        viewer.data().point_size = 10.0; // default is 30
        viewer.data().set_face_based(true);
    }

    return false;
}

int main(int argc, char *argv[])
{
  // Inline mesh of a cube
  /*V= (Eigen::MatrixXd(8,3)<<
    0.0,0.0,0.0,
    0.0,0.0,1.0,
    0.0,1.0,0.0,
    0.0,1.0,1.0,
    1.0,0.0,0.0,
    1.0,0.0,1.0,
    1.0,1.0,0.0,
    1.0,1.0,1.0).finished();
  F = (Eigen::MatrixXi(12,3)<<
    0,6,4,
    0,2,6,
    0,3,2,
    0,1,3,
    2,7,6,
    2,3,7,
    4,6,7,
    4,7,5,
    0,4,5,
    0,5,1,
    1,5,7,
    1,7,3).finished();*/

    // Load a surface mesh
    igl::readOBJ("../assets/cube.obj", V, F);
    //igl::readPLY("../assets/Armadillo.ply", V, F);
    
    //string filePath = "../assets/bunny.stl";
    //ifstream stlAscii(filePath);
    //if (stlAscii.is_open()) {
    //    igl::readSTL(stlAscii, V, F, N);
    //}

  // Tetrahedralize the interior
  igl::copyleft::tetgen::tetrahedralize(V, F, "pq1.414Y", TV, TT, TF);

  // Compute barycenters
  igl::barycenter(TV, TT, B);
  std::cout << "#row: " << B.rows() << "#col: " << B.cols() << std::endl;

  // Init the viewer
  igl::opengl::glfw::Viewer viewer;

  /////////////////////////////////////////////////////////////////////////
  //                               GUI                                   //
  /////////////////////////////////////////////////////////////////////////
 
  // Attach a menu plugin
  igl::opengl::glfw::imgui::ImGuiPlugin plugin;
  viewer.plugins.push_back(&plugin);
  igl::opengl::glfw::imgui::ImGuiMenu menu;
  plugin.widgets.push_back(&menu);

  // Add content to the default menu window
  menu.callback_draw_viewer_menu = [&]()
      {
          // Draw parent menu content, there're also keyboard shortcuts for the parent menu operations
          //menu.draw_viewer_menu();

          // Add new group
          if (ImGui::CollapsingHeader("Fracture Configuration Options:", ImGuiTreeNodeFlags_DefaultOpen))
          {
              // Expose an enumeration type
              enum Pattern { Uniform = 0, Radial, Brick, Grid, Slice };
              static Pattern type = Uniform;
              ImGui::Combo("Pattern", (int*)(&type), "Uniform\0Radial\0Brick\0Grid\0Slice\0\0");

              // We can also use a std::vector<std::string> defined dynamically
              static int num_choices = 3;
              static int idx_choice = 0;
              if (ImGui::InputInt("Num Nodes (for Node Placement)", &num_choices))
              {
                  num_choices = std::max(1, std::min(26, num_choices));
              }

              // Expose variable directly ...
              double doubleVariable = 0.1f;
              ImGui::InputDouble("Explode Amount", &doubleVariable, 0, 0, "%.4f");

              // ... or using a custom callback
              static bool partialFracture = false;
              if (ImGui::Checkbox("Partial Fracture", &partialFracture))
              {
                  // do something
                  std::cout << "partialFracture: " << std::boolalpha << partialFracture << std::endl;
              }

              // Add a button
              if (ImGui::Button("Pre-process VCAD", ImVec2(-1, 0)))
              {
                  std::cout << "Hello VACD\n";
              }
          }
      };

  // Plot the mesh
#if DEBUG
  viewer.callback_key_down = &key_down;
  key_down(viewer, '5', 0); // start off with 50% percentage of model cut
#else
  viewer.data().set_mesh(TV, TF);
  viewer.data().set_face_based(true);
  viewer.data().add_points(B, Eigen::RowVector3d(1, 0, 0));
  viewer.data().point_size = 10.0; // default is 30
#endif
  viewer.launch();
}
