#include <iostream>
#include <vector>
#include <random>

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

#include <chrono>

#include "tests.h"


using namespace std;

#define DEBUG       0
#define DEBUG_VORO  0
#define COCAD       0

// GUI
char gCurrKey = '0';
Eigen::IOFormat CleanFmt(4, 0, ", ", "\n", "[", "]");

enum MeshOp { Default = 0, Tet, Clip, Weld, Island, OBJ, Pipe};
static MeshOp gTestMode = Default;
double gExplodeAmt = 0.1f;

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
    viewer.data().add_points(convertToMatrixXd(gPoints), Eigen::RowVector3d(0, 1, 0));

    // Plot the voronoi cell boundary points as blue points
    Eigen::MatrixXd VV = convertToMatrixXd2D(gCellVertices);
#if DEBUG_VORO
    cout << "total cell verts: " << VV.rows() << endl;
    string sep = "\n----------------------------------------\n";
    Eigen::IOFormat cleanFmt(4, 0, ", ", "\n", "[", "]");
    cout << VV.format(cleanFmt) << sep << endl;
#endif
    viewer.data().add_points(VV, Eigen::RowVector3d(0, 0, 1));

    // Plot the voronoi cell edges as blue lines
    for (int cellIndex = 0; cellIndex < gCellEdges.size(); ++cellIndex) {
        const auto& edges = gCellEdges[cellIndex];
        const auto& vertices = gCellVertices[cellIndex];

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
    const auto& verts = gCellVertices[1];
    const auto& faces = gCellFaces[1];
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

// Core script to create splitted meshes 
// Needed data: V,F matrices of the whole mesh 
//              cellVertices, cellFaces, cellEdges from Voro Decomp
// Return:      std::vector<MeshConvex> a list of splitted convex mesh
void createMeshConvexs(std::vector<MeshConvex>& clippedMeshConvex, bool isMeshPrep = false) {
    clippedMeshConvex.clear();

    std::vector<Eigen::Vector3d> points;
    convertToVertArray(V, points);
    
    std::vector<std::vector<int>> faces;
    convertToFaceVertArray(F, faces);
    
    Surface_mesh tester_mesh;
    buildSMfromVF(points, faces, tester_mesh);
    
    MeshConvex tester{ points, faces, tester_mesh };
    calculateCentroid(tester, Eigen::Vector3d(0, 0, 0));
    
    vector<vector<Eigen::Vector3d>> cellVertices = gCellVertices;   // Vertices of each Voronoi cell
    vector<vector<vector<int>>> cellFaces = gCellFaces;
    
    // use an 8-cube pattern to decompose mesh first into compound
    if (isMeshPrep) {
        cellVertices.clear();
        cellFaces.clear();

        UnitCube cube;
        for (size_t i = 0; i < cube.direction.size(); ++i) {
            auto p = cube.points;
            for (auto& eachp : p) {
                eachp += Eigen::Vector3d(-0.5, -0.5, -0.5);
                eachp *= 0.5;
                eachp += cube.direction[i];
            }
            cellVertices.push_back(p);
            cellFaces.push_back(cube.faces);
        }

        //decomposeAABB(minCorner, maxCorner, cellVertices, cellFaces); //TODO: incorrect cell size, a cell covers the whole bunny
    }

    Pattern pattern(cellVertices, cellFaces, gCellEdges, minCorner, maxCorner); // we're not using gCellEdges rn, temp allowing to pass in edge info not matching verts&faces
    pattern.createCellsfromVoro();
    auto cells = pattern.getCells();
    
    for (auto const& c : cells) {
        spConvex result(new MeshConvex);
        clipConvexAgainstCell(tester, *c, result);
        calculateCentroid(*result, tester.centroid);
        clippedMeshConvex.push_back(*result);
    }
}

// Functional calls for moving each splitted mesh away 
// from the center of mass of the whole mesh by distance
// Needed data:     List of splitted mesh 
//                  distance of how far to move, larger, spreader
// Return:          MeshConvex stored vertices and faces of every splitted mesh  
MeshConvex splitMeshes(const std::vector<MeshConvex>& meshes, double distance) {
    std::vector<Eigen::Vector3d> final_vertices;
    std::vector<std::vector<int>> final_faces;
   
    for (auto& mesh : meshes) {
        int previous_verts = final_vertices.size();
        MeshConvex mesh_copy = mesh;
        translateMesh(mesh_copy, mesh_copy.centroid, distance);
        for (size_t i = 0; i < mesh_copy.faces.size(); i++) {
            mesh_copy.faces[i][0] += previous_verts;
            mesh_copy.faces[i][1] += previous_verts;
            mesh_copy.faces[i][2] += previous_verts;
        }
        final_vertices.insert(final_vertices.end(), mesh_copy.vertices.begin(), mesh_copy.vertices.end());
        final_faces.insert(final_faces.end(), mesh_copy.faces.begin(), mesh_copy.faces.end());
    }
    return MeshConvex{ final_vertices, final_faces };
}
MeshConvex splitMeshesPt(const std::vector<spConvex>& meshes, double distance) {
    std::vector<Eigen::Vector3d> final_vertices;
    std::vector<std::vector<int>> final_faces;

    for (auto& mesh : meshes) {
        int previous_verts = final_vertices.size();
        MeshConvex mesh_copy = *mesh;
        translateMesh(mesh_copy, mesh_copy.centroid, distance);
        for (size_t i = 0; i < mesh_copy.faces.size(); i++) {
            mesh_copy.faces[i][0] += previous_verts;
            mesh_copy.faces[i][1] += previous_verts;
            mesh_copy.faces[i][2] += previous_verts;
        }
        final_vertices.insert(final_vertices.end(), mesh_copy.vertices.begin(), mesh_copy.vertices.end());
        final_faces.insert(final_faces.end(), mesh_copy.faces.begin(), mesh_copy.faces.end());
    }
    return MeshConvex{ final_vertices, final_faces };
}

// Helper func to set viewer data
void setViewerMeshData(igl::opengl::glfw::Viewer& viewer,
    const std::vector<Eigen::Vector3d>& vertices,
    const std::vector<std::vector<int>>& faces) {
    auto V_temp = convertToMatrixXd(vertices);
    auto F_temp = convertToMatrixXi(faces);
    viewer.data().clear();
    viewer.data().set_mesh(V_temp, F_temp);
    viewer.data().set_face_based(true);
}

// Helper func, called every time a keyboard button is pressed
// Slice model at various percentage to view internal tetrahedral structure
bool key_down_tet(igl::opengl::glfw::Viewer& viewer, unsigned char key, int modifier) {
    using namespace Eigen;

    gCurrKey = key; // keep a global record

    if (key == '0')
    {
        viewer.data().clear();
    }
    if (key == '-')
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

    //drawDebugVisuals(viewer);

    return false;
}

// Helper func, called every time a keyboard button is pressed
// Display clipped convex in different indexed voronoi cell
bool key_down_clip(igl::opengl::glfw::Viewer& viewer, unsigned char key, int modifier) {
    gCurrKey = key; // keep a global record

    if (key >= '0' && key <= '9')
    {
        int cellIndex = int(key - '0');
        double scale = cellIndex * 0.01;
        auto mesh = splitMeshes(gClippedMeshConvex, scale);
        setViewerMeshData(viewer, mesh.vertices, mesh.faces);
    }

    drawDebugVisuals(viewer);

    return false;
}

// Helper func, called every time a keyboard button is pressed
// Draw the first convex within different indexed voronoi cell
bool key_down_weld(igl::opengl::glfw::Viewer& viewer, unsigned char key, int modifier) {
    using namespace Eigen;

    gCurrKey = key; // keep a global record

    if (key >= '0' && key <= '9')
    {
        int cellIndex = int(key - '0');
        auto mesh = gCells[cellIndex]->convexes[0]; // draw the first convex within this cell
        //auto mesh = gClippedMeshConvex[cellIndex]; // for testing cube positioning
        setViewerMeshData(viewer, mesh->vertices, mesh->faces);
    }

    drawDebugVisuals(viewer);

    return false;
}

// Helper func, called every time a keyboard button is pressed
// Input key='0~9': By default draws the first convex within a cell's compound
//      Input key = 'n/N': draws the other island(s) within the cell, if exists
bool key_down_island(igl::opengl::glfw::Viewer& viewer, unsigned char key, int modifier) {
    gCurrKey = key; // keep a global record
    if (key >= '0' && key <= '9')
    {
        int cellIndex = int(gCurrKey - '0');
        if (cellIndex >= gCompounds.size()) {
            std::cout << "There are " << gCompounds.size() << " compounds in total." << std::endl;
            viewer.data().clear();
            drawDebugVisuals(viewer);
            return false;
        }
        Compound temp = gCompounds[cellIndex];
        gCurrCompounds = islandDetection(temp);
        gCurrConvex = 0;
    }
    else if ((key == 'n' || key == 'N') && gCurrCompounds.size() > 0) {
        gCurrConvex = (gCurrConvex + 1) % gCurrCompounds.size();
    }

    if (gCurrCompounds.size() > 0) {
        Compound com = gCurrCompounds[gCurrConvex];
        std::vector<Eigen::Vector3d> final_vertices;
        std::vector<std::vector<int>> final_faces;
        // draw every convexes in current compound
        for (auto const& c : com.convexes) {
            int previous_verts = final_vertices.size();
            int previous_faces = final_faces.size();
            final_vertices.insert(final_vertices.end(), c->vertices.begin(), c->vertices.end());
            final_faces.insert(final_faces.end(), c->faces.begin(), c->faces.end());
            for (size_t i = previous_faces; i < final_faces.size(); i++) {
                final_faces[i][0] += previous_verts;
                final_faces[i][1] += previous_verts;
                final_faces[i][2] += previous_verts;
            }
            
        }
        setViewerMeshData(viewer, final_vertices, final_faces);
        drawDebugVisuals(viewer);
    }
    else {
        viewer.data().clear();
        drawDebugVisuals(viewer);
    }
    
    return false;
}

// Helper func, called every time a keyboard button is pressed
// Input key='0~9': By default draws the ith convex within a mesh compound
//      Input key = '-: draws all the convex hulls in a compound
bool key_down_obj(igl::opengl::glfw::Viewer& viewer, unsigned char key, int modifier) {
    using namespace Eigen;

    gCurrKey = key; // keep a global record
    if (key >= '0' && key <= '9')
    {
        int cellIndex = int(key - '0');
        setViewerMeshData(viewer, 
            ginitialConvexes.convexes[cellIndex]->vertices, 
            ginitialConvexes.convexes[cellIndex]->faces);
    }
    if (key == '-')
    {
        auto mesh = splitMeshesPt(ginitialConvexes.convexes, gExplodeAmt);
        setViewerMeshData(viewer, mesh.vertices, mesh.faces);
    }
    drawDebugVisuals(viewer);

    return false;
}

// Helper func, called every time a keyboard button is pressed
// Input key='0~9': By default draws the first convex within a cell's compound
//      Input key = 'n/N': draws the other island(s) within the cell, if exists
bool key_down_pipe(igl::opengl::glfw::Viewer& viewer, unsigned char key, int modifier) {
    using namespace Eigen;

    gCurrKey = key; // keep a global record
    std::vector<Eigen::Vector3d> final_vertices;
    std::vector<std::vector<int>> final_faces;

     //TODO: fix drawing all fracture pieces
    if (key == '-') {
        if (gFracturedConvex.size() < 1) {
            for (auto const& com : gCompounds) {
                for (auto const& c : com.convexes) {
                    gFracturedConvex.push_back(c);
                }
            }
        }

        auto mesh = splitMeshesPt(gFracturedConvex, gExplodeAmt);
        setViewerMeshData(viewer, mesh.vertices, mesh.faces);
        drawDebugVisuals(viewer);

        return false;
    }

    if (key >= '0' && key <= '9')
    {
        int cellIndex = int(key - '0');
        gCurrConvex = cellIndex; 
    }
    else if ((key == 'n' || key == 'N')) {
        gCurrConvex = gCurrConvex + 1;
    }

    if (gCurrConvex >= gCompounds.size()) {
        std::cout << "Error: gCurrConvex" << gCurrConvex << " exceed total compound number: " << gCompounds.size() << std::endl;
        return false;
    }

    // draw every convexes in current compound
    auto com = gCompounds[gCurrConvex];    
    for (auto const& c : com.convexes) {
        int previous_verts = final_vertices.size();
        int previous_faces = final_faces.size();
        final_vertices.insert(final_vertices.end(), c->vertices.begin(), c->vertices.end());
        final_faces.insert(final_faces.end(), c->faces.begin(), c->faces.end());
        for (size_t i = previous_faces; i < final_faces.size(); i++) {
            final_faces[i][0] += previous_verts;
            final_faces[i][1] += previous_verts;
            final_faces[i][2] += previous_verts;
        }

        setViewerMeshData(viewer, final_vertices, final_faces);
    }

    drawDebugVisuals(viewer);

    return false;
}

void switchTestMode(igl::opengl::glfw::Viewer& viewer)
{
    if (gTestMode == Default)
    {
        drawDebugVisuals(viewer);
    }
    if (gTestMode == Tet)
    {
        // Tetrahedralize the interior
        igl::copyleft::tetgen::tetrahedralize(V, F, "pq1.414Y", TV, TT, TF);
        // Compute barycenters
        igl::barycenter(TV, TT, B);

        viewer.callback_key_down = &key_down_tet;
        key_down_tet(viewer, '-', 0);   // '5', start off with 50% percentage of model cut;
                                        // '0', no model rendering
                                        // '-', full model
    }
    if (gTestMode == Clip)
    {
        createMeshConvexs(gClippedMeshConvex);

        viewer.callback_key_down = &key_down_clip;
        key_down_clip(viewer, '0', 0);
    }
    // Hardcoded 8 cubes as convex hulls
    if (gTestMode == Weld)
    {
        pre_test_welding(V, F);
        testWelding(gClippedMeshConvex, gCellVertices, gCellFaces, gCellEdges, gCells);

        viewer.callback_key_down = &key_down_weld;
        key_down_weld(viewer, '0', 0);
    }
    // Hardcoded 4 cubes (2 upper right row, 2 lower left row) as input compound
    // TODO: crashes in Release build
    if (gTestMode == Island)
    {
        testIsland(gClippedMeshConvex, gCellVertices, gCellFaces, gCellEdges, gCompounds);

        viewer.callback_key_down = &key_down_island;
        key_down_island(viewer, '0', 0);
    }
    if (gTestMode == OBJ)
    {
        testObj(gOBJPath, ginitialConvexes.convexes);
        
        viewer.callback_key_down = &key_down_obj;
        key_down_obj(viewer, '-', 0);
    }
    if (gTestMode == Pipe)
    {
        Pattern pattern(gCellVertices, gCellFaces, gCellEdges, minCorner, maxCorner);
        pattern.createCellsfromVoro();
        testPipeline(gOBJPath, pattern, gCompounds);
        
        viewer.callback_key_down = &key_down_pipe;
        key_down_pipe(viewer, '-', 0);
    }
}

int main(int argc, char *argv[])
{
    /////////////////////////////////////////////////////////////////////////
    //                         Load mesh                                   //
    /////////////////////////////////////////////////////////////////////////
    //string filePath = "../assets/obj/bunny.obj";//*"../assets/obj/cube.obj"; *//*".. /assets/bunny.off";*/ // "../assets/Armadillo.ply";
    string filePath = "../assets/obj/homer.obj";
    igl::readOBJ(filePath, V, F);
    //igl::readOFF(filePath, V, F);
     
    //igl::readPLY(filePath, V, F);
    //ifstream stlAscii(filePath);
    //if (stlAscii.is_open()) {
    //    igl::readSTL(stlAscii, V, F, N);
    //}

#if COCAD
  // run CoACD executable to decompose surface mesh into convex hulls
  LPCSTR applicationName = "..\\coacd.exe";
  char commandLine[] = "-i ..\\assets\\obj\\homer.obj -o ..\\assets\\results\\homer_out.obj -ro ..\\assets\\results\\remesh.obj";
  if (executeCommand(applicationName, commandLine))
  {
      cout << "CoACD executed successfully!" << endl;
  }
  else
  {
      cout << "CoACD execution failed!" << endl;
      return -1;
  }
#endif

  /////////////////////////////////////////////////////////////////////////
  //                         Voronoi decomposition                       //
  /////////////////////////////////////////////////////////////////////////
  // Get the mesh's bounding box
  minCorner = V.colwise().minCoeff();
  maxCorner = V.colwise().maxCoeff();
  std::cout << "Pattern min: \n" << minCorner << std::endl;
  std::cout << "Pattern max: \n" << maxCorner << std::endl;

  minCorner += Eigen::Vector3d(0, 0.3, 0);
  //minCorner -= Eigen::Vector3d(0.1, 0.1, 0.1);
  //maxCorner += Eigen::Vector3d(0.1, 0.1, 0.1);


#if VORO_LIB
  generateRandomPoints(gNumPoints, gPoints);
  convertEigenToVec(gPoints, gPointsVec);
  computeVoronoiCells(gPointsVec,
      vec3(minCorner.x(), minCorner.y(), minCorner.z()),
      vec3(maxCorner.x(), maxCorner.y(), maxCorner.z()),
      gCellVertices, gCellFaces, gCellEdges);
#endif

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
              if (ImGui::Combo("Test Mesh Operations", (int*)(&gTestMode), "Default\0Tet\0Clip\0Weld\0Island\0OBJ\0Pipe\0\0"))
              {
                  switchTestMode(viewer);
              }

              // We can also use a std::vector<std::string> defined dynamically
              if (ImGui::InputInt("Num Nodes (for Node Placement)", &gNumPoints))
              {
#if VORO_LIB
                  generateRandomPoints(gNumPoints, gPoints);
                  convertEigenToVec(gPoints, gPointsVec);
                  computeVoronoiCells(gPointsVec,
                      vec3(minCorner.x(), minCorner.y(), minCorner.z()),
                      vec3(maxCorner.x(), maxCorner.y(), maxCorner.z()),
                      gCellVertices, gCellFaces, gCellEdges);
#endif
                  key_down_tet(viewer, gCurrKey, 0);
              }

              // Expose variable directly ...
              if (ImGui::InputDouble("Explode Amount", &gExplodeAmt, 0, 0, "%.4f"))
              {
                  if (gTestMode == Clip && gClippedMeshConvex.size() > 0)
                  {
                      auto mesh = splitMeshes(gClippedMeshConvex, gExplodeAmt);
                      setViewerMeshData(viewer, mesh.vertices, mesh.faces);
                      drawDebugVisuals(viewer);
                  }
                  if (gTestMode == OBJ && ginitialConvexes.convexes.size() > 0)
                  {
                      auto mesh = splitMeshesPt(ginitialConvexes.convexes, gExplodeAmt);
                      setViewerMeshData(viewer, mesh.vertices, mesh.faces);
                      drawDebugVisuals(viewer);
                  }
                  if (gTestMode == Pipe && gFracturedConvex.size() > 0)
                  {
                      auto mesh = splitMeshesPt(gFracturedConvex, gExplodeAmt);
                      setViewerMeshData(viewer, mesh.vertices, mesh.faces);
                      drawDebugVisuals(viewer);
                  }
              }

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
  switchTestMode(viewer);

  Eigen::MatrixXd origin(1, 3); origin << 0, 0, 0;
  viewer.data().add_points(origin, Eigen::RowVector3d(1, 0, 0));

  viewer.launch();
}
