#ifndef MESH_CLIPPER
#define MESH_CLIPPER

#include <Eigen/Dense>
#include <Eigen/Sparse>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Plane_3.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/convex_hull_3.h>
#include <CGAL/Polygon_mesh_processing/triangulate_faces.h>
#include <CGAL/Polygon_mesh_processing/corefinement.h>
#include <CGAL/Polygon_mesh_processing/measure.h>
#include <CGAL/Polygon_mesh_processing/repair.h>
#include <CGAL/Polygon_mesh_processing/clip.h>
#include <CGAL/centroid.h>
#include <vector>
#include <memory>

typedef CGAL::Exact_predicates_inexact_constructions_kernel  K;
typedef CGAL::Polyhedron_3<K>                                Polyhedron_3;
typedef K::Point_3                                           Point_3;
typedef K::Vector_3                                          Vector_3;
typedef K::Plane_3                                           Plane_3;
typedef CGAL::Surface_mesh<Point_3>                          Surface_mesh;
typedef Surface_mesh::Vertex_index                                   vertex_descriptor;
typedef Surface_mesh::Face_index                                     face_descriptor;
namespace PMP = CGAL::Polygon_mesh_processing;

// representing a convex piece of the mesh, we only need vertices because the piece is convex
// and we can build the convex hull easily 
struct MeshConvex {
    std::vector<Eigen::Vector3d> vertices;  // Assuming convex shape is represented by its vertices
    std::vector<std::vector<int>> faces; 
    Surface_mesh convexMesh;                // convenient for cgal operations  
    Eigen::Vector3d centroid{0, 0, 0};      // center of mass of a set of points
};

// representing a convex piece of the fracture pattern
struct Cell {
    // could possibly add vertices and faces for future computations 
    std::vector<int> convexes; // Store convex index for intersected pieces 
    Surface_mesh cellMesh;            // convenient for cgal operations  
};

class Pattern 
{
    typedef std::vector<std::vector<Eigen::Vector3d>>  AllCellVertices;
    typedef std::vector<std::vector<std::vector<int>>> AllCellFaces;
    typedef std::vector<Eigen::MatrixXi>               AllCellEdges; 
    typedef std::shared_ptr<Cell>                      sPCell;

public:
    Pattern(AllCellVertices, AllCellFaces, AllCellEdges);
    ~Pattern() {};

    void createCellsfromVoro();
    std::vector<sPCell> getCells();

private:
    // original data outputed from Voronoi function of libigl
    AllCellVertices o_cellVertices;
    AllCellFaces o_cellFaces;
    AllCellEdges o_cellEdges;
    // data converted for mesh clipping and convex hull algorithm
    std::vector<sPCell> cells;

};

// convert vertices+faces data into triangulated Surface_mesh data
void buildSMfromVF(const std::vector<Eigen::Vector3d>&, const std::vector<std::vector<int>>&, Surface_mesh&);
// convert triangulated Surface_mesh data into vertices+faces data
void buildVFfromSM(const Surface_mesh&, std::vector<Eigen::Vector3d>&, std::vector<std::vector<int>>&);
// calculate the centroid of a MeshConvex
void calculateCentroid(MeshConvex& mesh, Eigen::Vector3d com); 
// move vertices in MeshConvex by scale in direction
void translateMesh(MeshConvex& mesh, Eigen::Vector3d direction, double scale);
MeshConvex clipConvexAgainstCell(const MeshConvex& convex, const Cell& cell);
#endif // !MESH_CLIPPER
