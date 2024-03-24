#include <iostream>
#include <vector>
#include <random>
#include "voro++.hh"
#include <igl/opengl/glfw/Viewer.h>

#include <igl/opengl/glfw/imgui/ImGuiPlugin.h>
#include <igl/opengl/glfw/imgui/ImGuiMenu.h>
#include <igl/opengl/glfw/imgui/ImGuiHelpers.h>

#include <igl/copyleft/tetgen/tetrahedralize.h>
#include <igl/readOBJ.h>
#include <igl/readPLY.h>
#include <igl/readSTL.h>
#include <igl/barycenter.h>

//#include <geogram/basic/common.h>
//#include <geogram/basic/logger.h>
//#include <geogram/basic/command_line.h>
//#include <geogram/mesh/mesh.h>
//#include <geogram/mesh/mesh_io.h>
//#include <geogram/voronoi/RVD.h>
//#include <geogram/delaunay/periodic_delaunay_3d.h>

//#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
//#include <CGAL/Delaunay_triangulation_3.h>
//#include <CGAL/Triangulation_vertex_base_with_info_3.h>
//#include <CGAL/point_generators_3.h>
//#include <CGAL/draw_triangulation_3.h>

//typedef CGAL::Exact_predicates_inexact_constructions_kernel         K;
//typedef CGAL::Triangulation_vertex_base_with_info_3<unsigned, K>    Vb;
//typedef CGAL::Triangulation_data_structure_3<Vb>                    Tds;
//typedef CGAL::Delaunay_triangulation_3<K, Tds>                      Delaunay;
//typedef Delaunay::Point                                             Point;
//typedef CGAL::Creator_uniform_3<double, Point>                      Creator;

using namespace std;

#define DEBUG 1

// Input polygon
Eigen::MatrixXd V; // #V by 3 matrix for vertices
Eigen::MatrixXi F; // matrix for face indices
Eigen::MatrixXd B; // matrix for barycenters
Eigen::MatrixXd N; // matrix for normals

// Tetrahedralized interior
Eigen::MatrixXd TV; // #TV by 3 matrix for vertex positions
Eigen::MatrixXi TT; // #TT by 4 matrix for tet face indices
Eigen::MatrixXi TF; // #TF by 3 matrix for triangle face indices ('f', else `boundary_facets` is called on TT)

// Voronoi diagram
bool periodic = false;
int numPoints = 3;

// Helper func, called every time a keyboard button is pressed
// Slice model at various percentage to view internal tetrahedral structure
bool key_down(igl::opengl::glfw::Viewer& viewer, unsigned char key, int modifier)
{
    using namespace Eigen;

    if (key == '0') return false;

    if (key >= '1' && key <= '9')
    {
        double t = double((key - '1') + 1) / 9.0;

        VectorXd v = B.col(2).array() - B.col(2).minCoeff();
        v /= v.col(0).maxCoeff();

        std::vector<int> s;

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

/**
 * \brief Initializes the Delaunay triangulation with
 *  random points.
 * \param[in] nb_points_in number of points
 * \param[delaunay] 
 */
// using namespace GEO;
//void init_random_points(int nb_points_in, PeriodicDelaunay3d &delaunay) {
//    index_t nb_points = index_t(nb_points_in);
//    nodes.resize(3 * nb_points);
//    for (index_t i = 0; i < 3 * nb_points; ++i) {
//        nodes[i] = Numeric::random_float64();
//    }
//    delaunay.set_vertices(nb_points, nodes.data());
//    delaunay.compute();
//}

void computeVoronoiCells(
    const vector<Eigen::Vector3d>& points,
    const Eigen::Vector3d& minCorner,
    const Eigen::Vector3d& maxCorner,
    vector<vector<Eigen::Vector3d>>& cellVertices, // TODO: user Eigen Matrix for 2d arrays
    vector<vector<int>>& cellFaces) {

    // Initialize container
    voro::container con(
        minCorner.x(), maxCorner.x(), 
        minCorner.y(), maxCorner.y(), 
        minCorner.z(), maxCorner.z(), 
        6, 6, 6,                // the number of grid blocks the container is divided into for computational efficiency.
        false, false, false,    // flags setting whether the container is periodic in x/y/z direction
        8);                     // allocate space for 8 particles in each computational block

    // Add particles to the container
    for (int i = 0; i < points.size(); i++) {
        con.put(i, points[i].x(), points[i].y(), points[i].z());
    }

    // Prepare cell computation
    voro::c_loop_all cl(con);
    voro::voronoicell_neighbor c;
    vector<double> fVerts; // Face vertices in global pos
    double x, y, z;

    if (cl.start()) do if (con.compute_cell(c, cl)) {
        // Extract cell vertices
        cl.pos(x, y, z);
        c.vertices(x, y, z, fVerts);
#ifdef DEBUG
        cout << "cell verts #: " << fVerts.size() << endl;
        for (int i = 0; i < fVerts.size(); i++)
        {
            cout << fVerts[i];
            if (i % 3 == 2)
            {
                cout << endl;
            }
            else {
                cout << ", ";
            }
        }
        cout << "\n" << endl;
#endif

        vector<Eigen::Vector3d> verts;
        // TOOD: do we want a triangulated verts array?
        for (int i = 0; i < fVerts.size(); i += 3) {
            verts.push_back(Eigen::Vector3d(fVerts[i], fVerts[i + 1], fVerts[i + 2]));
        }
        cellVertices.push_back(verts);

        // TODO: fix converting fVerts of doubles to faces of vertex indices
        // Extract cell faces (as vertex indices) for visualization
        vector<int> faces;
        vector<double>::iterator it = fVerts.begin();
        while (it != fVerts.end()) {
            int faceVerts = *(it++);
            faces.push_back(faceVerts); // Assuming libigl can handle these faces
            it += faceVerts * 3; // Skip the vertices of this face
        }
        cellFaces.push_back(faces);
    } while (cl.inc());
}

Eigen::MatrixXd convertToEigenMatrix(const vector<vector<Eigen::Vector3d>>& vectorOfVectors) {
    // First, calculate the total number of rows required in the matrix
    int totalRows = 0;
    for (const auto& vec : vectorOfVectors) {
        totalRows += vec.size();
    }

    // Create an Eigen::MatrixXd with the correct size
    Eigen::MatrixXd result(totalRows, 3); // 3 columns for x, y, z coordinates

    // Fill the matrix
    int currentRow = 0;
    for (const auto& vec : vectorOfVectors) {
        for (const auto& point : vec) {
            result(currentRow, 0) = point.x();
            result(currentRow, 1) = point.y();
            result(currentRow, 2) = point.z();
            currentRow++;
        }
    }

    return result;
}

Eigen::MatrixXd convertToMatrixXd(const std::vector<Eigen::Vector3d>& vec) {
    Eigen::MatrixXd mat(vec.size(), 3); // 3 columns for x, y, z
    for (size_t i = 0; i < vec.size(); ++i) {
        mat.row(i) = vec[i];
    }

    return mat;
}

int main(int argc, char *argv[])
{
    // Load a surface mesh
    string filePath = "../assets/cube.obj";// "../assets/bunny.stl"; // "../assets/Armadillo.ply"
    igl::readOBJ(filePath, V, F);
    //igl::readPLY(filePath, V, F);
    
    //ifstream stlAscii(filePath);
    //if (stlAscii.is_open()) {
    //    igl::readSTL(stlAscii, V, F, N);
    //}

  // Tetrahedralize the interior
  igl::copyleft::tetgen::tetrahedralize(V, F, "pq1.414Y", TV, TT, TF);

  // Compute barycenters
  igl::barycenter(TV, TT, B);
  std::cout << "barycenter #row: " << B.rows() << "#col: " << B.cols() << std::endl;

  // Voronoi decomposition
  std::vector<Eigen::Vector3d> points;

  // Get the mesh's bounding box
  Eigen::Vector3d minCorner = V.colwise().minCoeff();
  Eigen::Vector3d maxCorner = V.colwise().maxCoeff();

  // Generate random points
  // TODO: can look into https://doc.cgal.org/latest/Generator/Generator_2random_points_in_tetrahedral_mesh_3_8cpp-example.html#_a9
  std::random_device rd;  // Obtain random number from hardware and seed the generator
  std::mt19937 gen(rd());
  std::uniform_real_distribution<> disX(minCorner.x(), maxCorner.x());
  std::uniform_real_distribution<> disY(minCorner.y(), maxCorner.y());
  std::uniform_real_distribution<> disZ(minCorner.z(), maxCorner.z());

  for (int i = 0; i < numPoints; ++i) {
      double x = disX(gen);
      double y = disY(gen);
      double z = disZ(gen);
      points.push_back(Eigen::Vector3d(x, y, z)); 
  }

  //Delaunay dt(points.begin(), points.end());
  //cout << "dt.dimention = " << dt.dimension() << endl;

  vector<vector<Eigen::Vector3d>> cellVertices; // TODO: user Eigen Matrix for 2d arrays
  vector<vector<int>> cellFaces;
  computeVoronoiCells(points, minCorner, maxCorner, cellVertices, cellFaces);

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
  key_down(viewer, '0', 0); // start off with 50% percentage of model cut
#else
  viewer.data().set_mesh(TV, TF);
  viewer.data().set_face_based(true);
#endif
  // Corners of bounding box
  Eigen::MatrixXd V_box(8, 3);
  V_box <<
      minCorner(0), minCorner(1), minCorner(2),
      maxCorner(0), minCorner(1), minCorner(2),
      maxCorner(0), maxCorner(1), minCorner(2),
      minCorner(0), maxCorner(1), minCorner(2),
      minCorner(0), minCorner(1), maxCorner(2),
      maxCorner(0), minCorner(1), maxCorner(2),
      maxCorner(0), maxCorner(1), maxCorner(2),
      minCorner(0), maxCorner(1), maxCorner(2);

  // Edges of the bounding box
  Eigen::MatrixXi E_box(12, 2);
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
      7, 3;
  
  // Plot the corners of the bounding box as red points
  viewer.data().add_points(V_box, Eigen::RowVector3d(1, 0, 0));
  viewer.data().point_size = 30.0; // default is 30
  // Plot the edges of the bounding box
  for (unsigned i = 0; i < E_box.rows(); ++i)
  {
      viewer.data().add_edges
      (
          V_box.row(E_box(i, 0)),
          V_box.row(E_box(i, 1)),
          Eigen::RowVector3d(1, 0, 0)
      );
  }

  // Plot the voronoi cell nodes as green points
  viewer.data().add_points(convertToMatrixXd(points), Eigen::RowVector3d(0, 1, 0));
  // Plot the voronoi cell boundary points as blue points
  Eigen::MatrixXd VV(cellVertices[0].size() + cellVertices[1].size()+ cellVertices[2].size(), 3);
  for (int i = 0; i < cellVertices.size(); i++)
  {
      int start = VV.rows();
      cout << "# VV rows: " << start << endl;
      for (int j = 0; j < cellVertices[i].size(); j++)
      {
          Eigen::Vector3d pos = cellVertices[i][j];
          //VV(j + start, 0) = pos(0);
          //VV(j + start, 1) = pos(1);
          //VV(j + start, 2) = pos(2);
          cout << "VV row(" << j+start << ") = " << pos(0) << ", " << pos(1) << ", " << pos(2) << endl;
      }
  }
  viewer.data().point_size = 10.0;
  viewer.data().add_points(VV, Eigen::RowVector3d(0, 0, 1));
  viewer.launch();
}
