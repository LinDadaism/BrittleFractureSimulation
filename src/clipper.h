#pragma once

#include <Eigen/Dense>
#include <Eigen/Sparse>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Plane_3.h>
#include <CGAL/Bbox_3.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/convex_hull_3.h>
#include <CGAL/Polygon_mesh_processing/triangulate_faces.h>
#include <CGAL/Polygon_mesh_processing/corefinement.h>
#include <CGAL/Polygon_mesh_processing/measure.h>
#include <CGAL/Polygon_mesh_processing/repair.h>
#include <CGAL/Polygon_mesh_processing/clip.h>
#include <CGAL/Polygon_mesh_processing/bbox.h>
#include <CGAL/Polygon_mesh_processing/orientation.h>
#include <CGAL/centroid.h>
#include <CGAL/Rigid_triangle_mesh_collision_detection.h>
#include <CGAL/Convex_hull_3/dual/halfspace_intersection_3.h>
#include <CGAL/Convex_hull_3/dual/halfspace_intersection_with_constructions_3.h>
#include <CGAL/Polytope_distance_d.h>
#include <CGAL/Polytope_distance_d_traits_3.h>
#include <vector>
#include <thread>
#include <queue>
#include <mutex>
#include <condition_variable>
#include <algorithm>
#include <unordered_set>
#include <unordered_map>


#define VORO_LIB 1
#if VORO_LIB
    #include "../../../BFX/voro/include/voro++.hh"
#endif

typedef CGAL::Exact_predicates_inexact_constructions_kernel  K;
typedef CGAL::Polyhedron_3<K>                                Polyhedron_3;
typedef K::Point_3                                           Point_3;
typedef K::Vector_3                                          Vector_3;
typedef K::Plane_3                                           Plane_3;
typedef CGAL::Surface_mesh<Point_3>                          Surface_mesh;
typedef Surface_mesh::Vertex_index                           Vertex_descriptor;
typedef Surface_mesh::Face_index                             Face_descriptor;
#include <CGAL/Gmpzf.h>
typedef CGAL::Gmpzf ET;
typedef CGAL::Polytope_distance_d_traits_3<K, ET, double>    Traits;
typedef CGAL::Polytope_distance_d<Traits>                    Polytope_distance;
typedef CGAL::Rigid_triangle_mesh_collision_detection<Surface_mesh> ABtree;
namespace PMP = CGAL::Polygon_mesh_processing;
extern std::mutex patternMutex;             // global mutex for simplicity

// representing a convex piece of the mesh, we only need vertices because the piece is convex
// and we can build the convex hull easily 
struct MeshConvex {
    std::vector<Eigen::Vector3d> vertices;  // Assuming convex shape is represented by its vertices
    std::vector<std::vector<int>> faces; 
    Surface_mesh convexMesh;                // convenient for cgal operations  
    Eigen::Vector3d centroid{0, 0, 0};      // center of mass of a set of points
    double volume; 
};
typedef std::shared_ptr<MeshConvex>                          spConvex;

// representing a convex piece of the fracture pattern
struct Cell {
    int id;                                 // useful to track back all vertices and faces in Pattern 
    std::vector<spConvex> convexes;         // Store mesh convexes for intersected pieces 
    std::vector<Eigen::Vector3d> vertices;
    std::vector<std::vector<int>> faces;
    Surface_mesh cellMesh;                  // convenient for cgal operations  
};

struct Compound {
    std::vector<spConvex> convexes;         // a compound is consisted of bunch of convex pieces
    //Eigen::Vector3d minCorner{ 0, 0, 0 };   // min & max corners of compound's AABB
    //Eigen::Vector3d maxCorner{ 0, 0, 0 };
    Eigen::Vector3d centroid{ 0, 0, 0 };    // center of mass of a set of convex pieces
};
typedef std::shared_ptr<Compound>                           spCompound;

// Voronoi pattern typedefs
typedef std::vector<std::vector<Eigen::Vector3d>>  AllCellVertices;
typedef std::vector<std::vector<std::vector<int>>> AllCellFaces;
typedef std::vector<Eigen::MatrixXi>               AllCellEdges;
typedef std::pair<int, int>                        Edge;               // Edge represented by pair of vertex indices
typedef std::shared_ptr<Cell>                      spCell;

class Pattern
{
public:
    Pattern();
    Pattern(AllCellVertices, AllCellFaces, AllCellEdges, Eigen::Vector3d, Eigen::Vector3d);
    ~Pattern() {};

    void createCellsfromVoro();
    std::vector<spCell> getCells() const;
    AllCellVertices getVertices() const;
    AllCellFaces getFaces() const;
    Eigen::Vector3d getMin() const;
    Eigen::Vector3d getMax() const;
    int numCells() const;

    void setVertices(const AllCellVertices& verts);
    void setFaces(const AllCellFaces& faces);
    void setMin(const Eigen::Vector3d m); 
    void setMax(const Eigen::Vector3d m);

private:
    // original data outputed from Voronoi function of libigl
    AllCellVertices o_cellVertices;
    AllCellFaces o_cellFaces;
    AllCellEdges o_cellEdges;
    // data converted for mesh clipping and convex hull algorithm
    std::vector<spCell> cells;
    Eigen::Vector3d minCorner; 
    Eigen::Vector3d maxCorner;
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
// calculate the centroid of a Compound
Eigen::Vector3d calculateCentroidCompound(const std::vector<spConvex>& comp);
// move vertices in MeshConvex by scale in direction
void translateMesh(MeshConvex& mesh, Eigen::Vector3d direction, double scale);

// use a single cell to mesh clip a single convex piece
bool clipConvexAgainstCell(const MeshConvex& convex, const Cell& cell, spConvex& out_convex);
// Accelerate the clipping process using AABB tree 
void clipAABB(Compound& compound, Pattern& pattern);

// weld process for each cell in the pattern
void weldforPattern(Pattern& pattern);

// calculate the bounding box of any compound
void calculateBBox(const Compound& compound, Eigen::Vector3d& minCorner, Eigen::Vector3d& maxCorner);

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

// pre-fracture for the partial fracture
void preFracture(const Compound& compound, const Pattern& pattern, Compound& inside, Compound& outside);

#if VORO_LIB
void computeVoronoiCells(
    const std::vector<Eigen::Vector3d>& points,
    Eigen::Vector3d minCorner,
    Eigen::Vector3d maxCorner,
    AllCellVertices& cellVertices,  // TODO: use Eigen Matrix
    AllCellFaces& cellFaces,        // Each cell's faces by vertex indices
    AllCellEdges& cellEdges         // Each cell's edges
);
#endif