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

#include "clipper.h"

typedef std::pair<int, int> Edge; // Edge represented by pair of vertex indices

using namespace std;

#define DEBUG 1

// Input polygon
Eigen::MatrixXd V; // #V by 3 matrix for vertices
Eigen::MatrixXi F; // matrix for face indices
Eigen::MatrixXd B; // matrix for barycenters
Eigen::MatrixXd N; // matrix for normals
Eigen::Vector3d minCorner, maxCorner; // min and max corners of mesh's bounding box

// Tetrahedralized interior
Eigen::MatrixXd TV; // #TV by 3 matrix for vertex positions
Eigen::MatrixXi TT; // #TT by 4 matrix for tet face indices
Eigen::MatrixXi TF; // #TF by 3 matrix for triangle face indices ('f', else `boundary_facets` is called on TT)

// Voronoi diagram
bool periodic = false;
int numPoints = 8;
vector<Eigen::Vector3d> points;
vector<vector<Eigen::Vector3d>> cellVertices; // TODO: user Eigen Matrix for 2d arrays
vector<vector<vector<int>>> cellFaces;
vector<Eigen::MatrixXi> cellEdges;

// other global variables
char currKey = '0';

Eigen::MatrixXd convertToMatrixXd2D(const vector<vector<Eigen::Vector3d>>& vectorOfVectors) {
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
    for (int i = 0; i < vec.size(); ++i) {
        mat.row(i) = vec[i];
    }

    return mat;
}

void drawDebugVisuals(igl::opengl::glfw::Viewer& viewer) {
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
    viewer.data().point_size = 10.0; // default is 30, one size applies to all points visualization
    //viewer.data().add_points(V_box, Eigen::RowVector3d(1, 0, 0));
    
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
    Eigen::MatrixXd VV = convertToMatrixXd2D(cellVertices);
#ifdef DEBUG1
    cout << "total cell verts: " << VV.rows() << endl;
    string sep = "\n----------------------------------------\n";
    Eigen::IOFormat cleanFmt(4, 0, ", ", "\n", "[", "]");
    cout << VV.format(cleanFmt) << sep << endl;
#endif
    viewer.data().add_points(VV, Eigen::RowVector3d(0, 0, 1));

    // Plot the voronoi cell edges as blue lines
    for (int cellIndex = 0; cellIndex < cellEdges.size(); ++cellIndex) {
        const auto& edges = cellEdges[cellIndex];
        const auto& vertices = cellVertices[cellIndex];

        Eigen::MatrixXd V = convertToMatrixXd(vertices);
        Eigen::MatrixXd P1(edges.rows(), 3);
        Eigen::MatrixXd P2(edges.rows(), 3);

        for (int i = 0; i < edges.rows(); ++i) {
            int v1 = edges(i, 0);
            int v2 = edges(i, 1);

            P1.row(i) = V.row(v1);
            P2.row(i) = V.row(v2);
        }

        // Add edges to the viewer for this cell
        viewer.data().add_edges(P1, P2, Eigen::RowVector3d(0, 0, 1));
    }
}

// Helper func, called every time a keyboard button is pressed
// Slice model at various percentage to view internal tetrahedral structure
bool key_down(igl::opengl::glfw::Viewer& viewer, unsigned char key, int modifier) {
    using namespace Eigen;

    currKey = key; // keep a global record

    if (key == '0')
    {
        viewer.data().clear();
    }
    if (key == '*')
    {
        viewer.data().clear();
        viewer.data().set_mesh(V, F);
        viewer.data().set_face_based(true);
    }

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

    drawDebugVisuals(viewer);

    return false;
}

void generateRandomPoints(int numPoints, std::vector<Eigen::Vector3d>& points)
{
    points.clear();

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
}

void computeVoronoiCells(
    const vector<Eigen::Vector3d>& points,
    const Eigen::Vector3d& minCorner,
    const Eigen::Vector3d& maxCorner,
    vector<vector<Eigen::Vector3d>>& cellVertices, // TODO: use Eigen Matrix
    vector<vector<vector<int>>>& cellFaces, // Each cell's faces by vertex indices
    vector<Eigen::MatrixXi>& cellEdges // Each cell's edges
) {
    // fully clear nested vectors
    cellVertices.clear();
    cellFaces.clear();
    cellEdges.clear();

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
    double x, y, z;

    if (cl.start()) do if (con.compute_cell(c, cl)) {
        vector<double> cVerts;
        vector<int> neigh, fVerts;

        cl.pos(x, y, z);
        c.vertices(x, y, z, cVerts); // returns a vector of triplets (x,y,z)
        c.face_vertices(fVerts); // returns information about which vertices comprise each face:
                                 // It is a vector of integers with a specific format: 
                                 // the first entry is a number k corresponding to the number of vertices
                                 // making up a face, and this is followed k additional entries 
                                 // describing which vertices make up this face. 
                                 // e.g. (3, 16, 20, 13) would correspond to a triangular face linking
                                 // vertices 16, 20, and 13 together.
        c.neighbors(neigh); // returns the neighboring particle IDs corresponding to each face

#ifdef DEBUG1
        cout << "cell verts #: " << cVerts.size() / 3 << endl;
        for (int i = 0; i < cVerts.size(); i++)
        {
            cout << cVerts[i];
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
        // Process vertices
        vector<Eigen::Vector3d> verts;
        for (int i = 0; i < cVerts.size(); i += 3) {
            verts.push_back(Eigen::Vector3d(cVerts[i], cVerts[i + 1], cVerts[i + 2]));
        }

        // Process cell faces and edges
        vector<vector<int>> faces; // faces<face vertex indices>
        vector<Edge> edges;
        for (int i = 0; i < fVerts.size(); i += fVerts[i] + 1) {
            int n = fVerts[i]; // Number of vertices for this face
            vector<int> face;

            for (int k = 1; k <= n; ++k) {
                face.push_back(fVerts[i + k]);

                // Extract edges
                int currVertex = fVerts[i + k];
                int nextVertex = fVerts[i + ((k % n) + 1)];
                Edge edge(std::min(currVertex, nextVertex), std::max(currVertex, nextVertex));

                if (std::find(edges.begin(), edges.end(), edge) == edges.end()) {
                    edges.push_back(edge); // Add unique edge
                }
            }

            faces.push_back(face);
        }

        cellVertices.push_back(verts);
        cellFaces.push_back(faces);

        // TODO: move to a helper func
        // Convert edgesVector to Eigen::MatrixXi for current cell
        Eigen::MatrixXi edgesMatrix(edges.size(), 2);
        for (int i = 0; i < edges.size(); ++i) {
            edgesMatrix(i, 0) = edges[i].first;
            edgesMatrix(i, 1) = edges[i].second;
        }
        cellEdges.push_back(edgesMatrix);

    } while (cl.inc());
}
void testFunc() {
    typedef CGAL::Exact_predicates_inexact_constructions_kernel  K;
    typedef CGAL::Polyhedron_3<K>                     Polyhedron_3;
    typedef K::Point_3                                Point_3;
    typedef CGAL::Surface_mesh<Point_3>               Surface_mesh;
    std::vector<Point_3> points;
    Point_3 p;
    points.push_back(Point_3(0, 0, 0));
    points.push_back(Point_3(1, 0, 0));
    points.push_back(Point_3(0, 1, 0));
    points.push_back(Point_3(1, 1, 0));
    points.push_back(Point_3(0, 0, 1));
    points.push_back(Point_3(1, 0, 1));
    points.push_back(Point_3(0, 1, 1));
    points.push_back(Point_3(1, 1, 1));
    // define polyhedron to hold convex hull
    Polyhedron_3 poly;
    // compute convex hull of non-collinear points
    CGAL::convex_hull_3(points.begin(), points.end(), poly);
    std::cout << "The convex hull contains " << poly.size_of_vertices() << " vertices" << std::endl;
    Surface_mesh sm;
    CGAL::convex_hull_3(points.begin(), points.end(), sm);
    std::cout << "The convex hull contains " << num_vertices(sm) << " vertices" << std::endl;
}

int main(int argc, char *argv[])
{
    /////////////////////////////////////////////////////////////////////////
    //                         Load mesh                                   //
    /////////////////////////////////////////////////////////////////////////
    string filePath = "../assets/cube.obj";// "../assets/bunny.stl"; // "../assets/Armadillo.ply"
    igl::readOBJ(filePath, V, F);
    //igl::readPLY(filePath, V, F);
    
    //ifstream stlAscii(filePath);
    //if (stlAscii.is_open()) {
    //    igl::readSTL(stlAscii, V, F, N);
    //}
    testFunc();
  // Tetrahedralize the interior
  igl::copyleft::tetgen::tetrahedralize(V, F, "pq1.414Y", TV, TT, TF);

  // Compute barycenters
  igl::barycenter(TV, TT, B);
  std::cout << "barycenter #row: " << B.rows() << "#col: " << B.cols() << std::endl;

    /////////////////////////////////////////////////////////////////////////
    //                         Voronoi decomposition                       //
    /////////////////////////////////////////////////////////////////////////

  // Get the mesh's bounding box
  minCorner = V.colwise().minCoeff();
  maxCorner = V.colwise().maxCoeff();

  generateRandomPoints(numPoints, points);
  computeVoronoiCells(points, minCorner, maxCorner, cellVertices, cellFaces, cellEdges);

  /*std::vector<Eigen::Vector3d> verts; 
  verts.push_back(Eigen::Vector3d(1., 0., 0.));
  verts.push_back(Eigen::Vector3d(0., 1., 0.));
  verts.push_back(Eigen::Vector3d(0., 0., 1.));*/
  /*MeshConvex meshconvex{ verts };
  buildConvexHull(meshconvex);*/
  /////////////////////////////////////////////////////////////////////////
  //                               GUI                                   //
  /////////////////////////////////////////////////////////////////////////
  // Init the viewer
  igl::opengl::glfw::Viewer viewer;

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
              if (ImGui::InputInt("Num Nodes (for Node Placement)", &numPoints))
              {
                  generateRandomPoints(numPoints, points);
                  computeVoronoiCells(points, minCorner, maxCorner, cellVertices, cellFaces, cellEdges);
                  key_down(viewer, currKey, 0);
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

  /////////////////////////////////////////////////////////////////////////
  //                        Visualization                                //
  /////////////////////////////////////////////////////////////////////////
  // Plot the tet mesh
#if DEBUG
  viewer.callback_key_down = &key_down;
  key_down(viewer, '5', 0); // '5', start off with 50% percentage of model cut;
                            // '0', no model rendering
                            // '*', full model
#else
  viewer.data().set_mesh(TV, TF);
  viewer.data().set_face_based(true);
#endif

  viewer.launch();
}
