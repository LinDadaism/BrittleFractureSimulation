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
#include <algorithm>
#include <unordered_set>
#include <unordered_map>

typedef CGAL::Exact_predicates_inexact_constructions_kernel  K;
typedef CGAL::Polyhedron_3<K>                                Polyhedron_3;
typedef K::Point_3                                           Point_3;
typedef K::Vector_3                                          Vector_3;
typedef K::Plane_3                                           Plane_3;
typedef CGAL::Surface_mesh<Point_3>                          Surface_mesh;
typedef Surface_mesh::Vertex_index                           Vertex_descriptor;
typedef Surface_mesh::Face_index                             Face_descriptor;
namespace PMP = CGAL::Polygon_mesh_processing;

// representing a convex piece of the mesh, we only need vertices because the piece is convex
// and we can build the convex hull easily 
struct MeshConvex {
    std::vector<Eigen::Vector3d> vertices;  // Assuming convex shape is represented by its vertices
    std::vector<std::vector<int>> faces; 
    Surface_mesh convexMesh;                // convenient for cgal operations  
    Eigen::Vector3d centroid{0, 0, 0};      // center of mass of a set of points
    double volume; 
    std::unordered_set<int> group;               // data structure used in island detection
};
typedef std::shared_ptr<MeshConvex>                          spConvex;

// representing a convex piece of the fracture pattern
struct Cell {
    int id; // useful to track back all vertices and faces in Pattern 
    std::vector<spConvex> convexes; // Store mesh convexes for intersected pieces 
    Surface_mesh cellMesh;            // convenient for cgal operations  
};

struct Compound {
    std::vector<spConvex> convexes; // a compound is consisted of bunch of convex pieces 
};
typedef std::shared_ptr<Compound>                           spCompound;

class Pattern 
{
public:
    typedef std::vector<std::vector<Eigen::Vector3d>>  AllCellVertices;
    typedef std::vector<std::vector<std::vector<int>>> AllCellFaces;
    typedef std::vector<Eigen::MatrixXi>               AllCellEdges;
    typedef std::shared_ptr<Cell>                      sPCell;

    Pattern(AllCellVertices, AllCellFaces, AllCellEdges);
    ~Pattern() {};

    void createCellsfromVoro();
    std::vector<sPCell> getCells();
    AllCellVertices getVertices(); 
    AllCellFaces getFaces(); 

private:
    // original data outputed from Voronoi function of libigl
    AllCellVertices o_cellVertices;
    AllCellFaces o_cellFaces;
    AllCellEdges o_cellEdges;
    // data converted for mesh clipping and convex hull algorithm
    std::vector<sPCell> cells;
    // possible data to be added 
    // 1. the center of the pattern or the bounding box of Pattern 
    // 2. transformation matrices to move the pattern around
    // 3. scal

};

// convert vertices+faces data into triangulated Surface_mesh data
void buildSMfromVF(const std::vector<Eigen::Vector3d>&, const std::vector<std::vector<int>>&, Surface_mesh&);
// convert triangulated Surface_mesh data into vertices+faces data
void buildVFfromSM(const Surface_mesh&, std::vector<Eigen::Vector3d>&, std::vector<std::vector<int>>&);
// calculate the centroid of a MeshConvex
void calculateCentroid(MeshConvex& mesh, Eigen::Vector3d com); 
// move vertices in MeshConvex by scale in direction
void translateMesh(MeshConvex& mesh, Eigen::Vector3d direction, double scale);
// use a single cell to mesh clip a single convex piece
bool clipConvexAgainstCell(const MeshConvex& convex, const Cell& cell, spConvex& out_convex);
// weld process for each cell in the pattern
void weldforPattern(Pattern& pattern);
// create a plane from 3 points
Eigen::Vector4d createPlane(Eigen::Vector3d p1, Eigen::Vector3d p2, Eigen::Vector3d p3);
// the boring bit - injecting a hash specialisation into the std:: namespace
// but let's derive from boost's hash class, which is much better
// in that it allows easy hashing using free functions
namespace std {
    template<> struct hash<::Eigen::Vector4d> : boost::hash<::Eigen::Vector4d> {};
}
// island detection algorithm 
std::vector<Compound> islandDetection(Compound& old_compound);
// pipeline of fracture algorithm 
std::vector<Compound> fracturePipeline(Compound& compound, Pattern& pattern); 
#endif // !MESH_CLIPPER
