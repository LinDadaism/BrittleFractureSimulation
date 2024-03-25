#include "clipper.h"

MeshConvex clipConvexAgainstPlane(const MeshConvex& convex, const Plane& plane) {
    std::vector<Eigen::Vector3d> newVertices;

    for (size_t i = 0; i < convex.vertices.size(); ++i) {
        size_t next = (i + 1) % convex.vertices.size();
        Eigen::Vector3d currentVertex = convex.vertices[i];
        Eigen::Vector3d nextVertex = convex.vertices[next];

        bool currentInside = !plane.isPointOnPositiveSide(currentVertex);
        bool nextInside = !plane.isPointOnPositiveSide(nextVertex);

        if (currentInside) {
            newVertices.push_back(currentVertex);
        }

        if (currentInside != nextInside) {
            float t = (plane.distance - plane.normal.dot(currentVertex)) / plane.normal.dot(nextVertex - currentVertex);
            Eigen::Vector3d intersection = currentVertex + t * (nextVertex - currentVertex);
            newVertices.push_back(intersection);
        }
    }

    return MeshConvex{ newVertices };
}

MeshConvex clipConvexAgainstCell(const MeshConvex& convex, const Cell& cell) {
    MeshConvex curMesh;
    curMesh.vertices = convex.vertices;
    for (const auto& plane : cell.planes) {
        curMesh = clipConvexAgainstPlane(curMesh, plane);
    }
    return curMesh;
}

//std::pair<Eigen::MatrixXd, Eigen::MatrixXi> buildConvexHull(const MeshConvex& convex) {
//    typedef CGAL::Exact_predicates_inexact_constructions_kernel  K;
//    typedef CGAL::Polyhedron_3<K>                     Polyhedron_3;
//    typedef K::Point_3                                Point_3;
//    typedef CGAL::Surface_mesh<Point_3>               Surface_mesh;
//
//    std::vector<Point_3> points;
//    for (const auto& v : convex.vertices) {
//        points.push_back(Point_3(v[0], v[1], v[2]));
//    }
//    // define polyhedron to hold convex hull
//    Polyhedron_3 poly;
//    // compute convex hull of non-collinear points
//    CGAL::convex_hull_3(points.begin(), points.end(), poly);
//    std::cout << "The convex hull contains " << poly.size_of_vertices() << " vertices" << std::endl;
//    Surface_mesh sm;
//    CGAL::convex_hull_3(points.begin(), points.end(), sm);
//    std::cout << "The convex hull contains " << num_vertices(sm) << " vertices" << std::endl;
//    Eigen::MatrixXd a; 
//    Eigen::MatrixXi b; 
//    return std::make_pair(a, b);
//}

