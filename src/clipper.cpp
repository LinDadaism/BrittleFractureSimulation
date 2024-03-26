#include "clipper.h"

Pattern::Pattern(AllCellVertices v, AllCellFaces f, AllCellEdges e):o_cellVertices(v), o_cellFaces(f), o_cellEdges(e)
{}

std::vector<Cell> Pattern::getCells() {
    return cells;
}

void Pattern::createCellsfromVoro() {
    for (size_t i = 0; i < o_cellFaces.size(); ++i) {
        auto curCell = o_cellFaces[i];
        std::vector<Plane> planes; 
        for (size_t j = 0; j < curCell.size(); ++j) {
            int n1 = curCell[j][0];
            int n2 = curCell[j][1];
            int n3 = curCell[j][2];
            Eigen::Vector3d p1 = o_cellVertices[i][n1]; 
            Eigen::Vector3d p2 = o_cellVertices[i][n2];
            Eigen::Vector3d p3 = o_cellVertices[i][n3];
            Eigen::Vector3d v1 = p2 - p1; 
            Eigen::Vector3d v2 = p3 - p1; 
            Eigen::Vector3d normal = v1.cross(v2); 
            double d = normal.dot(p1); 
            planes.push_back(Plane{ -normal, -d});
        }
        cells.push_back(Cell{ planes, std::vector<MeshConvex>{} });
    }
}

MeshConvex clipConvexAgainstPlane(const MeshConvex& convex, const Plane& plane) {
    std::vector<Eigen::Vector3d> newVertices;

    for (const auto& f : convex.faces) {
        for (size_t i = 0; i < f.size(); ++i) {
            int cur = f[i]; 
            int next = f[(i + 1) % f.size()]; 
            Eigen::Vector3d currentVertice = convex.vertices[cur]; 
            Eigen::Vector3d nextVertice = convex.vertices[next];
            bool currentInside = !plane.isPointOnPositiveSide(currentVertice);
            bool nextInside = !plane.isPointOnPositiveSide(nextVertice);
            if (currentInside) {
                newVertices.push_back(currentVertice);
            }
            if (currentInside != nextInside) {
                float t = (plane.distance - plane.normal.dot(currentVertice)) / plane.normal.dot(nextVertice - currentVertice);
                Eigen::Vector3d intersection = currentVertice + t * (nextVertice - currentVertice);
                newVertices.push_back(intersection);
            }
        }
    }
    // create the convex hull
    std::vector<Point_3> points;
    for (const auto& v : newVertices) {
        points.push_back(Point_3(v.x(), v.y(), v.z()));
    }
    Surface_mesh sm; 
    CGAL::convex_hull_3(points.begin(), points.end(), sm);
    std::vector<Eigen::Vector3d> resultVertices;
    std::vector<std::vector<int>> resultFaces; 
    for (const auto& v : sm.vertices()) {
        const Point_3& p = sm.point(v);
        resultVertices.push_back(Eigen::Vector3d(p.x(), p.y(), p.z()));
    }
    for (const auto& f : sm.faces()) {
        std::vector<int> face;
        for (auto vd : vertices_around_face(sm.halfedge(f), sm)) {
            face.push_back(vd.idx());
        }
        resultFaces.push_back(face);
    }
    return MeshConvex{ resultVertices, resultFaces };
}

MeshConvex clipConvexAgainstCell(const MeshConvex& convex, const Cell& cell) {
    MeshConvex curMesh = convex;
    for (const auto& plane : cell.planes) {
        curMesh = clipConvexAgainstPlane(curMesh, plane);
    }
    return curMesh;
}
