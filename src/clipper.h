#ifndef MESH_CLIPPER
#define MESH_CLIPPER

#include <Eigen/Dense>
#include <Eigen/Sparse>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/convex_hull_3.h>
#include <vector>
#include <memory>

typedef CGAL::Exact_predicates_inexact_constructions_kernel  K;
typedef CGAL::Polyhedron_3<K>                     Polyhedron_3;
typedef K::Point_3                                Point_3;
typedef K::Vector_3                               Vector_3;
typedef CGAL::Surface_mesh<Point_3>               Surface_mesh;

struct Plane {
    Eigen::Vector3d normal;
    double distance;  // Distance from the origin

    // Assuming a point is on the positive side if it's in the direction of the normal
    bool isPointOnPositiveSide(const Eigen::Vector3d& point) const {
        return normal.dot(point) - distance > 0;
    }
};

// representing a convex piece of the mesh, we only need vertices because the piece is convex
// and we can build the convex hull easily 
struct MeshConvex {
    std::vector<Eigen::Vector3d> vertices;  // Assuming convex shape is represented by its vertices
    std::vector<std::vector<int>> faces; 
};

// representing a convex piece of the fracture pattern
struct Cell {
    std::vector<Plane> planes;
    std::vector<MeshConvex> convexes; // Store intersected pieces 
};

class Pattern 
{
    typedef std::vector<std::vector<Eigen::Vector3d>>  AllCellVertices;
    typedef std::vector<std::vector<std::vector<int>>> AllCellFaces;
    typedef std::vector<Eigen::MatrixXi>               AllCellEdges; 

public:
    Pattern(AllCellVertices, AllCellFaces, AllCellEdges);
    ~Pattern() {};

    void createCellsfromVoro();
    std::vector<Cell> getCells(); 

private:
    // original data outputed from Voronoi function of libigl
    AllCellVertices o_cellVertices;
    AllCellFaces o_cellFaces;
    AllCellEdges o_cellEdges;
    // data converted for mesh clipping and convex hull algorithm
    std::vector<Cell> cells;

};
MeshConvex clipConvexAgainstPlane(const MeshConvex& convex, const Plane& plane);
MeshConvex clipConvexAgainstCell(const MeshConvex& convex, const Cell& cell);

#endif // !MESH_CLIPPER
