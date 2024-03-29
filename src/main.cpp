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
#include <igl/readOFF.h>
#include <igl/barycenter.h>

#include "clipper.h"
//#include <igl/copyleft/cgal/convex_hull.h>

typedef std::pair<int, int> Edge; // Edge represented by pair of vertex indices

using namespace std;

#define DEBUG 0
#define DEBUG_VORO 0

// Input polygon
Eigen::MatrixXd V; // #V by 3 matrix for vertices
Eigen::MatrixXi F; // matrix for face indices
Eigen::MatrixXd B; // matrix for barycenters
Eigen::MatrixXd N; // matrix for normals
Eigen::Vector3d minCorner, maxCorner; // min and max corners of mesh's bounding box

Eigen::MatrixXd meshV;
Eigen::MatrixXi meshF;

// Tetrahedralized interior
Eigen::MatrixXd TV; // #TV by 3 matrix for vertex positions
Eigen::MatrixXi TT; // #TT by 4 matrix for tet face indices
Eigen::MatrixXi TF; // #TF by 3 matrix for triangle face indices ('f', else `boundary_facets` is called on TT)

// Voronoi diagram
bool periodic = false;
int numPoints = 10;
vector<Eigen::Vector3d> points;

// TODO: encapsulate class Compound
vector<vector<Eigen::Vector3d>> cellVertices;   // Vertices of each Voronoi cell
vector<vector<vector<int>>> cellFaces;          // Faces of each Voronoi cell represented by vertex indices
vector<Eigen::MatrixXi> cellEdges;              // Edges of each Voronoi cell represented by vertex indices
K::Point_3 centerOfMass; 

// other global variables
char currKey = '0';
Eigen::IOFormat CleanFmt(4, 0, ", ", "\n", "[", "]");

std::vector<MeshConvex> clippedMeshConvex;

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

Eigen::MatrixXi convertToMatrixXi(const std::vector<std::vector<int>>& vec) {
    Eigen::MatrixXi mat(vec.size(), 3); // 3 columns for x, y, z
    for (int i = 0; i < vec.size(); ++i) {
        mat(i, 0) = vec[i][0];
        mat(i, 1) = vec[i][1];
        mat(i, 2) = vec[i][2];
    }

    return mat;
}

void convertToFaceVertArray(const Eigen::MatrixXi& faceVertsMat, vector<vector<int>>& faceVerts)
{
    for (int i = 0; i < faceVertsMat.rows(); ++i) {
        // Temporary vector to hold the vertex indices for the current face
        std::vector<int> indices;

        for (int j = 0; j < faceVertsMat.cols(); ++j) {
            indices.push_back(faceVertsMat(i, j));
        }
        faceVerts.push_back(indices);
    }
}

void convertToVertArray(const Eigen::MatrixXd& vertMatrix, vector<Eigen::Vector3d>& vertVector)
{
    for (int i = 0; i < vertMatrix.rows(); ++i) {
        Eigen::Vector3d vertex = vertMatrix.row(i);
        vertVector.push_back(vertex);
    }
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
#if DEBUG_VORO
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

#if DEBUG_VORO
    // Plot the voronoi cell faces
    // Assuming we're visualizing the ith cell
    const auto& verts = cellVertices[1];
    const auto& faces = cellFaces[1];
    for (auto& f : faces)
    {
        cout << "# fverts: " << f.size() << endl;
    }
    // Convert vertices to Eigen::MatrixXd
    Eigen::MatrixXd CV(verts.size(), 3);
    for (size_t i = 0; i < verts.size(); ++i) {
        CV.row(i) = verts[i];
    }

    // Convert faces to Eigen::MatrixXi
    // Note: This assumes all faces are triangles. If not, you'll need to triangulate faces.
    size_t totalFaces = 0;
    for (const auto& face : faces) {
        totalFaces += face.size() - 2; // Simple triangulation for non-triangular faces
    }
    Eigen::MatrixXi CF(totalFaces, 3);
    size_t currentRow = 0;
    for (const auto& face : faces) {
        // Assuming the face is a simple polygon that can be triangulated by a fan.
        for (size_t i = 1; i + 1 < face.size(); ++i) {
            CF.row(currentRow++) << face[0], face[i], face[i + 1];
        }
    }
    cout << "CF:\n" << CF.format(CleanFmt) << endl;
    viewer.data().set_mesh(CV, CF);
    viewer.data().set_face_based(true);
#endif
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
            std::reverse(face.begin(), face.end());
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

MeshConvex testFunc2(int cellIndex) {
    /*std::vector<Eigen::Vector3d> points = { Eigen::Vector3d(0.0,0.0,0.0),
                                            Eigen::Vector3d(0.0,0.0,1.0),
                                            Eigen::Vector3d(0.0,1.0,0.0),
                                            Eigen::Vector3d(0.0,1.0,1.0),
                                            Eigen::Vector3d(1.0,0.0,0.0),
                                            Eigen::Vector3d(1.0,0.0,1.0),
                                            Eigen::Vector3d(1.0,1.0,0.0),
                                            Eigen::Vector3d(1.0,1.0,1.0) };
    for (auto& p : points) {
        p += Eigen::Vector3d(-0.5, -0.5, -0.5);
        p *= 1;
    }
    std::vector<std::vector<int>> faces = { {0,6,4},{0,2,6},{0,3,2},{0,1,3},
                                            {2,7,6},{2,3,7},{4,6,7},{4,7,5},
                                            {0,4,5},{0,5,1},{1,5,7},{1,7,3} };*/
    std::vector<Eigen::Vector3d> points;
    convertToVertArray(V, points);
    std::vector<std::vector<int>> faces; 
    convertToFaceVertArray(F, faces);
    Surface_mesh tester_mesh; 
    buildSMfromVF(points, faces, tester_mesh);
    MeshConvex tester{ points, faces, tester_mesh};
    Pattern pattern(cellVertices, cellFaces, cellEdges);
    pattern.createCellsfromVoro();
    Cell cell = *pattern.getCells()[cellIndex];
    auto result = clipConvexAgainstCell(tester, cell);
    calculateCentroid(result, Eigen::Vector3d(0, 0, 0));
    translateMesh(result, result.centroid, 0.0);
    return result;
}

// Draw every cell 
MeshConvex testFunc3() {
    // Convert format from V,F matrices
    std::vector<Eigen::Vector3d> points;
    convertToVertArray(V, points);
    std::vector<std::vector<int>> faces;
    convertToFaceVertArray(F, faces);
    // Build SM of the whole mesh
    Surface_mesh tester_mesh;
    buildSMfromVF(points, faces, tester_mesh);
    MeshConvex tester{ points, faces, tester_mesh };
    // Calculate the centroid of the whole mesh use for visualization
    calculateCentroid(tester, Eigen::Vector3d(0, 0, 0));
    // Construct Pattern from Voronoi Decomposition
    Pattern pattern(cellVertices, cellFaces, cellEdges);
    pattern.createCellsfromVoro();
    auto cells = pattern.getCells();
    std::vector<Eigen::Vector3d> final_vertices; 
    std::vector<std::vector<int>> final_faces; 
    // Calculate the center of mass for each clipped mesh and 
    // move it away from the center of mass of the whole mesh
    for (auto const& c : cells) {
        int previous_verts = final_vertices.size();
        auto result = clipConvexAgainstCell(tester, *c);
        calculateCentroid(result, tester.centroid);
        translateMesh(result, result.centroid, 0.03);
        for (size_t i = 0; i < result.faces.size(); i++) {
            result.faces[i][0] += previous_verts;
            result.faces[i][1] += previous_verts;
            result.faces[i][2] += previous_verts;
        }
        final_vertices.insert(final_vertices.end(), result.vertices.begin(), result.vertices.end());
        final_faces.insert(final_faces.end(), result.faces.begin(), result.faces.end());
    }

    return MeshConvex{final_vertices, final_faces};
}

// Core script to create splitted meshes 
// Needed data: V,F matrices of the whole mesh 
//              cellVertices, cellFaces, cellEdges from Voro Decomp
// Return:      std::vector<MeshConvex> a list of splitted convex mesh
std::vector<MeshConvex> createMeshes() {
    std::vector<Eigen::Vector3d> points;
    convertToVertArray(V, points);
    std::vector<std::vector<int>> faces;
    convertToFaceVertArray(F, faces);
    Surface_mesh tester_mesh;
    buildSMfromVF(points, faces, tester_mesh);
    MeshConvex tester{ points, faces, tester_mesh };
    calculateCentroid(tester, Eigen::Vector3d(0, 0, 0));
    Pattern pattern(cellVertices, cellFaces, cellEdges);
    pattern.createCellsfromVoro();
    auto cells = pattern.getCells();
    std::vector<MeshConvex> results; 
    for (auto const& c : cells) {
        auto result = clipConvexAgainstCell(tester, *c);
        calculateCentroid(result, tester.centroid);
        results.push_back(result);
    }
    return results;
}

// Functional calls for moving each splitted mesh away 
// from the center of mass of the whole mesh by distance
// Needed data:     List of splitted mesh 
//                  distance of how far to move, larger, spreader
// Return:          MeshConvex stored vertices and faces of every splitted mesh  
MeshConvex splitMeshes(const std::vector<MeshConvex> meshes, double distance) {
    std::vector<Eigen::Vector3d> final_vertices;
    std::vector<std::vector<int>> final_faces;
    for (auto& mesh : meshes) {
        int previous_verts = final_vertices.size();
        MeshConvex mesh_copy = mesh;
        translateMesh(mesh_copy, mesh_copy.centroid, distance);
        for (size_t i = 0; i < mesh.faces.size(); i++) {
            mesh_copy.faces[i][0] += previous_verts;
            mesh_copy.faces[i][1] += previous_verts;
            mesh_copy.faces[i][2] += previous_verts;
        }
        final_vertices.insert(final_vertices.end(), mesh_copy.vertices.begin(), mesh_copy.vertices.end());
        final_faces.insert(final_faces.end(), mesh_copy.faces.begin(), mesh_copy.faces.end());
    }
    return MeshConvex{ final_vertices, final_faces };
}

// Helper func, called every time a keyboard button is pressed
// Slice model at various percentage to view internal tetrahedral structure
bool key_down_clip(igl::opengl::glfw::Viewer& viewer, unsigned char key, int modifier) {
    using namespace Eigen;

    currKey = key; // keep a global record

    if (key >= '0' && key <= '9')
    {
        int cellIndex = int(key - '0');
        //auto mesh = testFunc2(cellIndex);
        //auto mesh = testFunc3();
        double scale = cellIndex * 0.01;
        auto mesh = splitMeshes(clippedMeshConvex, scale);
        auto V_temp = convertToMatrixXd(mesh.vertices);
        auto F_temp = convertToMatrixXi(mesh.faces);
        viewer.data().clear();
        viewer.data().set_mesh(V_temp, F_temp);
        viewer.data().set_face_based(true);
    }

    drawDebugVisuals(viewer);

    return false;
}

int main(int argc, char *argv[])
{
    /////////////////////////////////////////////////////////////////////////
    //                         Load mesh                                   //
    /////////////////////////////////////////////////////////////////////////
    string filePath = /*"../assets/cube.obj";*/ "../assets/bunny.off"; // "../assets/Armadillo.ply"
    //igl::readOBJ(filePath, V, F);
    igl::readOFF(filePath, V, F);
     
    //igl::readPLY(filePath, V, F);
    //ifstream stlAscii(filePath);
    //if (stlAscii.is_open()) {
    //    igl::readSTL(stlAscii, V, F, N);
    //}
    //igl::copyleft::cgal::convex_hull(meshV, V, F);

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

  clippedMeshConvex = createMeshes();
  
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
  viewer.data().clear();
  viewer.data().set_mesh(V, F);
  viewer.data().set_face_based(true);
  drawDebugVisuals(viewer);
#if DEBUG
  viewer.callback_key_down = &key_down;
  key_down(viewer, '0', 0); // '5', start off with 50% percentage of model cut;
                            // '0', no model rendering
                            // '*', full model
#else
  viewer.callback_key_down = &key_down_clip;
  key_down_clip(viewer, '0', 0);
  Eigen::MatrixXd origin(1, 3); origin << 0, 0, 0;
  viewer.data().add_points(origin, Eigen::RowVector3d(1, 0, 0));
#endif

  viewer.launch();
}
