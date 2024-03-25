#ifndef MESH_CLIPPER
#define MESH_CLIPPER

#include <Eigen/Dense>
#include <Eigen/Sparse>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/convex_hull_3.h>
#include <vector>
#include <fstream>


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
};

// representing a convex piece of the fracture pattern
struct Cell {
    std::vector<Plane> planes;
};

MeshConvex clipConvexAgainstPlane(const MeshConvex& convex, const Plane& plane);
MeshConvex clipConvexAgainstCell(const MeshConvex& convex, const Cell& cell);
//std::pair<Eigen::MatrixXd, Eigen::MatrixXi> buildConvexHull(const MeshConvex& convex);

//Mesh buildCGALMesh()
#endif // !MESH_CLIPPER
