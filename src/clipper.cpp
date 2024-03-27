#include "clipper.h"

Pattern::Pattern(AllCellVertices v, AllCellFaces f, AllCellEdges e):o_cellVertices(v), o_cellFaces(f), o_cellEdges(e)
{}

std::vector<Pattern::sPCell> Pattern::getCells() {
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
            planes.push_back(Plane{ normal, d});
        }
        Surface_mesh sm; 
        buildSMfromVF(o_cellVertices[i], o_cellFaces[i], sm);
        sPCell cellptr(new Cell{ planes, std::vector<int>{}, sm });
        cells.push_back(cellptr);
    }
}

void buildSMfromVF(const std::vector<Eigen::Vector3d>& vertices, const std::vector<std::vector<int>>& faces, Surface_mesh& sm) {
    sm.clear();
    // put in all vertices 
    for (const auto& v : vertices) {
        sm.add_vertex(Point_3(v.x(), v.y(), v.z()));
    }
    // put in all faces 
    for (const auto& f : faces) {
        std::vector<vertex_descriptor> indices; 
        for (int i : f) {
            indices.push_back(vertex_descriptor(i));
        }
        sm.add_face(indices);
    }
    // triangualte all non-triangle faces
    PMP::triangulate_faces(sm);
}

void buildVFfromSM(const Surface_mesh& sm, std::vector<Eigen::Vector3d>& vertices, std::vector<std::vector<int>>& faces) {
    vertices.clear(); 
    faces.clear();
    for (const auto& v : sm.vertices()) {
        const Point_3& p = sm.point(v);
        vertices.push_back(Eigen::Vector3d(p.x(), p.y(), p.z()));
    }
    for (const auto& f : sm.faces()) {
        std::vector<int> face;
        for (auto vd : vertices_around_face(sm.halfedge(f), sm)) {
            face.push_back(vd.idx());
        }
        faces.push_back(face);
    }
}

void calculateCentroid(MeshConvex& mesh, Eigen::Vector3d center) {
    auto sm = mesh.convexMesh;
    K::Vector_3 com(0.0, 0.0, 0.0); 
    K::Vector_3 old_center(center.x(), center.y(), center.z());
    int count = 0; 
    for (auto v : sm.vertices()) {
        com = com + (sm.point(v) - CGAL::ORIGIN - old_center);
        ++count;
    }
    if (count != 0) {
        com = com / count;
    }
    mesh.centroid = Eigen::Vector3d(com.x(), com.y(), com.z());
}

void translateMesh(MeshConvex& mesh, Eigen::Vector3d direction, double scale) {
    for (size_t i = 0; i < mesh.vertices.size(); i++) {
        mesh.vertices[i] += direction.normalized() * scale;
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

MeshConvex clipConvexAgainstPlane2(const MeshConvex& convex, const Plane& plane) {
    MeshConvex convex_copy = convex; 
    Plane_3 p(plane.normal.x(), plane.normal.y(), plane.normal.z(), plane.distance); 
    Surface_mesh convexm = convex_copy.convexMesh;
    bool check = PMP::clip(convexm, p, PMP::parameters::clip_volume(true).use_compact_clipper(true));
    convexm.collect_garbage();
    std::vector<Eigen::Vector3d> vertices; 
    std::vector<std::vector<int>> faces; 
    buildVFfromSM(convexm, vertices, faces);
    return MeshConvex{ vertices, faces, convexm };
}


MeshConvex clipConvexAgainstCell(const MeshConvex& convex, const Cell& cell) {
    MeshConvex curMesh = convex;
    for (const auto& plane : cell.planes) {
        curMesh = clipConvexAgainstPlane(curMesh, plane);
    }
    return curMesh;
}

MeshConvex clipConvexAgainstCell2(const MeshConvex& convex, const Cell& cell) {
    // use the copy to do the clipping as cgal clipping modifies on the original meshes
    MeshConvex convex_copy = convex;
    Cell       cell_copy = cell; 
    Surface_mesh convexm = convex_copy.convexMesh; 
    Surface_mesh cellm   = cell_copy.cellMesh;
    bool a = PMP::does_self_intersect(convexm);
    bool b = PMP::does_self_intersect(cellm);
    //while (a) {
    //    Surface_mesh::Halfedge_index hole_he_index;
    //    for (auto h : convexm.halfedges()) {
    //        if (convexm.is_border(h)) {
    //            hole_he_index = h;
    //            break; // Found a boundary halfedge, break the loop
    //        }
    //    }
    //    std::vector<Surface_mesh::Face_index>  patched_faces;
    //    PMP::triangulate_hole(mesh,
    //        hole_he_index,
    //        std::back_inserter(patched_faces),
    //        CGAL::parameters::vertex_point_map(mesh.points()).geom_traits(Kernel()));
    //}
    //PMP::clip(cellm, convexm, PMP::parameters::clip_volume(true), PMP::parameters::use_compact_clipper(true));
    PMP::clip(convexm, cellm, PMP::parameters::clip_volume(true).use_compact_clipper(true));
    //Surface_mesh result; 
    //PMP::corefine_and_compute_intersection(convexm, cellm, result);
    convexm.collect_garbage();
    std::vector<Eigen::Vector3d> vertices; 
    std::vector<std::vector<int>> faces; 
    buildVFfromSM(convexm, vertices, faces);
    return MeshConvex{ vertices, faces, convexm };
}