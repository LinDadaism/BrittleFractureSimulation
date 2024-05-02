#include "clipper.h"
#include <iostream>
#include <chrono>

#define MULTITHREAD 1
#if MULTITHREAD
    #define CGAL_HAS_THREADS
    #define BOOST_HAS_THREADS
#endif
#define INTERSECT_DISTANCE 1e-5
std::mutex patternMutex;

Pattern::Pattern() {}

Pattern::Pattern(AllCellVertices v, AllCellFaces f, AllCellEdges e, Eigen::Vector3d minc, Eigen::Vector3d maxc):
    o_cellVertices(v), o_cellFaces(f), o_cellEdges(e), minCorner(minc), maxCorner(maxc)
{}

std::vector<spCell> Pattern::getCells() const {
    return cells;
}

AllCellVertices Pattern::getVertices() const {
    return o_cellVertices;
}

AllCellFaces Pattern::getFaces() const {
    return o_cellFaces;
}

int Pattern::numCells() const {
    return cells.size();
}

Eigen::Vector3d Pattern::getMin() const {
    return minCorner;
}

Eigen::Vector3d Pattern::getMax() const {
    return maxCorner;
}

void Pattern::setVertices(const AllCellVertices& verts) {
    o_cellVertices.clear();

    for (int i = 0; i < verts.size(); i++) {
        auto& cell = verts[i];

        std::vector<Eigen::Vector3d> cellVerts;
        for (int j = 0; j < cell.size(); j++) {
            cellVerts.push_back(cell[j]);
        }
        o_cellVertices.push_back(cellVerts);
    }
}

void Pattern::setFaces(const AllCellFaces& faces) {
    o_cellFaces.clear();

    for (int i = 0; i < faces.size(); i++) {
        auto& cell = faces[i];

        std::vector<std::vector<int>> cellVertIds;
        for (int j = 0; j < cell.size(); j++) {
            auto& face = cell[j];

            std::vector<int> faceVertIds;
            for (int k = 0; k < face.size(); k++) {
                faceVertIds.push_back(face[k]);
            }
            cellVertIds.push_back(faceVertIds);
        }
        o_cellFaces.push_back(cellVertIds);
    }
}

void Pattern::setMin(const Eigen::Vector3d m) {
    minCorner = m; 
}
void Pattern::setMax(const Eigen::Vector3d m) {
    maxCorner = m;
}

void Pattern::createCellsfromVoro() {
    // for every cell create Cell and push to container
    for (int i = 0; i < o_cellFaces.size(); ++i) {
        auto curCell = o_cellFaces[i];
        Surface_mesh sm; 
        buildSMfromVF(o_cellVertices[i], o_cellFaces[i], sm);
        spCell cellptr(new Cell{ i, std::vector<spConvex>{}, o_cellVertices[i], o_cellFaces[i], sm });
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
        std::vector<Vertex_descriptor> indices; 
        for (int i : f) {
            indices.push_back(Vertex_descriptor(i));
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
    // old center represents the com of the whole mesh and we want 
    // to calculate the center of mass of a convex piece from the 
    // whole mesh
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

Eigen::Vector3d calculateCentroidCompound(const std::vector<spConvex>& comp) {
    Eigen::Vector3d totalWeightedCentroid(0, 0, 0);
    double totalVolume = 0;

    for (const auto& convex : comp) {
        totalWeightedCentroid += convex->centroid * convex->volume;
        totalVolume += convex->volume;
    }
    Eigen::Vector3d compoundCentroid = totalWeightedCentroid / totalVolume;

    return compoundCentroid;
}

void translateMesh(MeshConvex& mesh, Eigen::Vector3d direction, double scale) {
    auto d = direction.normalized();
    for (size_t i = 0; i < mesh.vertices.size(); i++) {
        mesh.vertices[i] += d * scale;
    }
}

void calculateBBox(const Compound& compound, Eigen::Vector3d& minCorner, Eigen::Vector3d& maxCorner) {
    if (compound.convexes.size() < 1) {
        std::cerr << "Compound is empty" << std::endl;
    }
    minCorner = Eigen::Vector3d(INFINITY, INFINITY, INFINITY);
    maxCorner = Eigen::Vector3d(-INFINITY, -INFINITY, -INFINITY);
    for (const auto& convex : compound.convexes) {
        auto box = PMP::bbox(convex->convexMesh);
        minCorner[0] = min(minCorner[0], box.xmin());
        minCorner[1] = min(minCorner[1], box.ymin());
        minCorner[2] = min(minCorner[2], box.zmin());
        maxCorner[0] = max(maxCorner[0], box.xmax());
        maxCorner[1] = max(maxCorner[1], box.ymax());
        maxCorner[2] = max(maxCorner[2], box.zmax());
    }
}

/*previous code of mesh clipping from scratch, temporarily saved 
* for possible future reference 
*/
//MeshConvex clipConvexAgainstPlane(const MeshConvex& convex, const Plane& plane) {
//    std::vector<Eigen::Vector3d> newVertices;
//
//    for (const auto& f : convex.faces) {
//        for (size_t i = 0; i < f.size(); ++i) {
//            int cur = f[i]; 
//            int next = f[(i + 1) % f.size()]; 
//            Eigen::Vector3d currentVertice = convex.vertices[cur]; 
//            Eigen::Vector3d nextVertice = convex.vertices[next];
//            bool currentInside = !plane.isPointOnPositiveSide(currentVertice);
//            bool nextInside = !plane.isPointOnPositiveSide(nextVertice);
//            if (currentInside) {
//                newVertices.push_back(currentVertice);
//            }
//            if (currentInside != nextInside) {
//                float t = (plane.distance - plane.normal.dot(currentVertice)) / plane.normal.dot(nextVertice - currentVertice);
//                Eigen::Vector3d intersection = currentVertice + t * (nextVertice - currentVertice);
//                newVertices.push_back(intersection);
//            }
//        }
//    }
//    // create the convex hull
//    std::vector<Point_3> points;
//    for (const auto& v : newVertices) {
//        points.push_back(Point_3(v.x(), v.y(), v.z()));
//    }
//    Surface_mesh sm; 
//    CGAL::convex_hull_3(points.begin(), points.end(), sm);
//    std::vector<Eigen::Vector3d> resultVertices;
//    std::vector<std::vector<int>> resultFaces; 
//    for (const auto& v : sm.vertices()) {
//        const Point_3& p = sm.point(v);
//        resultVertices.push_back(Eigen::Vector3d(p.x(), p.y(), p.z()));
//    }
//    for (const auto& f : sm.faces()) {
//        std::vector<int> face;
//        for (auto vd : vertices_around_face(sm.halfedge(f), sm)) {
//            face.push_back(vd.idx());
//        }
//        resultFaces.push_back(face);
//    }
//    return MeshConvex{ resultVertices, resultFaces };
//}

// Remeber to first detect collision then use this function.
// The clipped convex is guranteed to have volume greater than 0  
bool clipConvexAgainstCell(const MeshConvex& convex, const Cell& cell, spConvex& out_convex) {
    // use the copy to do the clipping as cgal clipping modifies on the original meshes
    // could do in place modification if that's faster
    MeshConvex convex_copy = convex;
    Cell       cell_copy = cell; 
    Surface_mesh convexm = convex_copy.convexMesh; 
    Surface_mesh cellm   = cell_copy.cellMesh;
  
    //Possible code to get rid of degenerated or self-intersect mesh 
    //commented for now since we only accept convex sound mesh in theory
    //bool a = PMP::does_self_intersect(convexm);
    //bool b = PMP::does_self_intersect(cellm);
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
    PMP::clip(convexm, cellm, PMP::parameters::clip_volume(true).use_compact_clipper(true));
    convexm.collect_garbage();
    double volume = PMP::volume(convexm);
    std::vector<Eigen::Vector3d> vertices;
    std::vector<std::vector<int>> faces;
    // double check that they actually intersected 
    if (volume > 0) {
        buildVFfromSM(convexm, vertices, faces);
    }
    else {
        return false;
    }
    spConvex result(new MeshConvex{ vertices, faces, convexm, {0, 0, 0}, volume });
    out_convex = result;
    return true;
}

#if MULTITHREAD
void workerClipAABB(int cellID, int cell_num, Compound& compound, Pattern& pattern, std::vector<std::pair<size_t, bool>> intersects) {
    Eigen::Vector3d compCentroid = compound.centroid;
    for (const auto& inter : intersects) {
        // for all intersections that is between a cell and a convex
        if (inter.first >= cell_num) {
            Surface_mesh convexm;
            CGAL::copy_face_graph(compound.convexes[inter.first - cell_num]->convexMesh, convexm);
            Surface_mesh cellm;
            CGAL::copy_face_graph(pattern.getCells()[cellID]->cellMesh, cellm);
            // clipping between mesh and mesh 
            //PMP::clip(convexm, cellm, PMP::parameters::clip_volume(true).use_compact_clipper(true));
            PMP::corefine_and_compute_intersection(convexm, cellm, convexm);
            convexm.collect_garbage();
            double volume = PMP::volume(convexm);
            std::vector<Eigen::Vector3d> vertices;
            std::vector<std::vector<int>> faces;
            // double check that they actually intersected 
            if (volume > 0) {
                buildVFfromSM(convexm, vertices, faces);
                spConvex result(new MeshConvex{ vertices, faces, convexm, {0, 0, 0}, volume });
                calculateCentroid(*result, Eigen::Vector3d(0, 0, 0));
                pattern.getCells()[cellID]->convexes.push_back(result);
            }
            
        }
    }
}
#endif

void clipAABB(Compound& compound, Pattern& pattern) {
    ABtree tree;
    int cell_num = pattern.getCells().size();
    int convex_num = compound.convexes.size();
    Eigen::Vector3d compCentroid = compound.centroid;

    // put all surface meshes into AABB tree 
    for (const auto& cell : pattern.getCells()) {
        // assume all meshes 1 connected pieces
        tree.add_mesh(cell->cellMesh, PMP::parameters::apply_per_connected_component(false));
    }
    for (const auto& convex : compound.convexes) {
        tree.add_mesh(convex->convexMesh, PMP::parameters::apply_per_connected_component(false));
    }
    std::vector<std::vector<std::pair<size_t, bool>>> intersection_list;
    for (int i = 0; i < cell_num; ++i) {
        std::vector<std::pair<size_t, bool>> intersects = tree.get_all_intersections_and_inclusions(i);
        intersection_list.push_back(intersects);
    }
    // iterate each cell and find all intersections 
#if MULTITHREAD
    std::vector<std::thread> threads;
    for (int i = 0; i < cell_num; ++i) {
        std::thread myThread(workerClipAABB, i, cell_num, compound, pattern, intersection_list[i]);
        threads.push_back(std::move(myThread));
    }
    for (auto& thread : threads) {
        thread.join();
    }
#else
    for (int i = 0; i < cell_num; ++i) {
        std::vector<std::pair<size_t, bool>> intersects = tree.get_all_intersections_and_inclusions(i);
        for (const auto& inter : intersects) {
            // for all intersections that is between a cell and a convex
            if (inter.first >= cell_num) {
                Surface_mesh convexm; 
                CGAL::copy_face_graph(compound.convexes[inter.first - cell_num]->convexMesh, convexm);
                Surface_mesh cellm; 
                CGAL::copy_face_graph(pattern.getCells()[i]->cellMesh, cellm);
                //clipping between mesh and mesh 
                auto start = std::chrono::high_resolution_clock::now();
                //PMP::clip(convexm, cellm, PMP::parameters::clip_volume(true).use_compact_clipper(true));
                PMP::corefine_and_compute_intersection(convexm, cellm, convexm);
                auto stop = std::chrono::high_resolution_clock::now();
                auto duration = std::chrono::duration_cast<std::chrono::microseconds>(stop - start);
                convexm.collect_garbage();
                double volume = PMP::volume(convexm);
                std::vector<Eigen::Vector3d> vertices;
                std::vector<std::vector<int>> faces;
                // double check that they actually intersected 
                if (volume > 0) {
                    buildVFfromSM(convexm, vertices, faces);
                }
                spConvex result(new MeshConvex{ vertices, faces, convexm, {0, 0, 0}, volume });
                calculateCentroid(*result, compCentroid);
                pattern.getCells()[i]->convexes.push_back(result);
            }
        }
    }
#endif
    // DEBUG
    //std::cout << "Clipping times after AABB tree: " << counter << std::endl;
    //std::cout << "Including times after AABB tree: " << includer << std::endl;
    //std::cout << "Total clipping time after AABB tree: " << total_time.count() << std::endl;

}


// make sure to run this function after clipping such that 
// cell.convexes contain clipped convex pieces 
void weldforPattern(Pattern& pattern) {
    for (auto& cell : pattern.getCells()) {
        auto cellsm = cell->cellMesh;
        double total_volume = PMP::volume(cellsm);
        double sum_volume = 0; 
        for (auto& convex : cell->convexes) {
            sum_volume += convex->volume;
        }
        if (sum_volume - total_volume >= 0) {
            // clear all convex pieces and replace with 1 piece
            // of the same shape and size of the cell
            cell->convexes.clear();
            auto id = cell->id; 
            std::vector<Eigen::Vector3d> vertices;
            std::vector<std::vector<int>> faces;
            buildVFfromSM(cellsm, vertices, faces);
            spConvex spconvex(new MeshConvex{ vertices, faces, cellsm, {0, 0, 0}, total_volume});
            cell->convexes.push_back(spconvex);
        }
    }
}

Eigen::Vector4d createPlane(Eigen::Vector3d p1, Eigen::Vector3d p2, Eigen::Vector3d p3) {
    Eigen::Vector3d v(p2.x() - p1.x(), p2.y() - p1.y(), p2.z() - p1.z());
    Eigen::Vector3d w(p3.x() - p1.x(), p3.y() - p1.y(), p3.z() - p1.z());
    auto A = (v.y() * w.z()) - (v.z() * w.y());
    auto B = (v.z() * w.x()) - (v.x() * w.z());
    auto C = (v.x() * w.y()) - (v.y() * w.x());
    auto D = -(A * p1.x() + B * p1.y() + C * p1.z());
    return Eigen::Vector4d(A, B, C, D);
}

//island detection using intersection 
std::vector<Compound> islandDetection(Compound& old_compound) {
    int N = old_compound.convexes.size();
    ABtree tree;
    for (auto& com : old_compound.convexes) tree.add_mesh(com->convexMesh, PMP::parameters::apply_per_connected_component(false));
    std::vector<Compound> results;
    std::unordered_set<int> cur_convs;
    std::unordered_set<int> final_convs;
    for (int l = 0; l < N; ++l) final_convs.insert(l);
    for (int i = 0; i < N; ++i) {
        if (cur_convs == final_convs) break;
        if (cur_convs.find(i) != cur_convs.end()) continue;
        std::unordered_set<int> searched{};
        std::queue<int> q;
        searched.insert(i);
        q.push(i);
        while (!q.empty()) {
            int current = q.front();
            q.pop();
            std::unordered_set<int> intersected;
            for (auto& intersection : tree.get_all_intersections_and_inclusions(current)) {
                intersected.insert(intersection.first);
            }
            for (auto& inter : intersected) {
                // not found
                if (searched.find(inter) == searched.end()) {
                    searched.insert(inter);
                    q.push(inter);
                }
            }
        }
        std::vector<std::shared_ptr<MeshConvex>> new_convexes;
        for (auto& s : searched) {
            new_convexes.push_back(old_compound.convexes[s]);
            cur_convs.insert(s);
        }
        results.push_back(Compound{ new_convexes });
    }
    
    return results;
}
//////////////////////////////////////////////////////////////////////////////
// The alternative island detection algorithm that use the distance between 
// each convexes which is much slower 
//std::vector<Compound> islandDetection(Compound& old_compound) {
//    int N = old_compound.convexes.size();
//    std::vector<Compound> results;
//    std::unordered_set<int> cur_convs;
//    std::unordered_set<int> final_convs;
//    for (int l = 0; l < N; ++l) final_convs.insert(l);
//    for (int i = 0; i < N; ++i) {
//        if (cur_convs == final_convs) break;
//        if (cur_convs.find(i) != cur_convs.end()) continue;
//        std::unordered_set<int> searched{};
//        std::queue<int> q;
//        searched.insert(i);
//        q.push(i);
//        while (!q.empty()) {
//            int current = q.front();
//            q.pop();
//            std::unordered_set<int> intersected;
//            for (int k = 0; k < N; ++k) {
//                if (current != k
//                    && searched.find(k) == searched.end()
//                    && cur_convs.find(k) == cur_convs.end()) {
//                    std::vector<Point_3> P;
//                    std::vector<Point_3> Q;
//                    auto sm_P = old_compound.convexes[current]->convexMesh;
//                    auto sm_Q = old_compound.convexes[k]->convexMesh;
//                    for (auto& v : sm_P.vertices()) P.push_back(sm_P.point(v));
//                    for (auto& v : sm_Q.vertices()) Q.push_back(sm_Q.point(v));
//                    Polytope_distance pd(P.begin(), P.end(), Q.begin(), Q.end());
//                    double sd = CGAL::to_double(pd.squared_distance_numerator()) /
//                        CGAL::to_double(pd.squared_distance_denominator()); 
//                    if (sd < INTERSECT_DISTANCE) intersected.insert(k);
//                }
//            }
//            for (auto& inter : intersected) {
//                // not found
//                if (searched.find(inter) == searched.end()) {
//                    searched.insert(inter);
//                    q.push(inter);
//                }
//            }
//        }
//        std::vector<std::shared_ptr<MeshConvex>> new_convexes;
//        for (auto& s : searched) {
//            new_convexes.push_back(old_compound.convexes[s]);
//            cur_convs.insert(s);
//        }
//        results.push_back(Compound{ new_convexes });
//    }
//
//    return results;
//}
////////////////////////////////////////////////////////////////////////////////

// pre-fracture for partial fracture
void preFracture(const Compound& compound, const Pattern& pattern, Compound& inside, Compound& outside) {
    ///////////////////////////////////////////////////
    // first create a cube sm for a given pattern
    Surface_mesh sm; 
    auto mi = pattern.getMin(); 
    auto ma = pattern.getMax();
    // put in all vertices 
    sm.add_vertex(Point_3(mi.x(), mi.y(), mi.z()));
    sm.add_vertex(Point_3(mi.x(), mi.y(), ma.z()));
    sm.add_vertex(Point_3(mi.x(), ma.y(), mi.z()));
    sm.add_vertex(Point_3(mi.x(), ma.y(), ma.z()));
    sm.add_vertex(Point_3(ma.x(), mi.y(), mi.z()));
    sm.add_vertex(Point_3(ma.x(), mi.y(), ma.z()));
    sm.add_vertex(Point_3(ma.x(), ma.y(), mi.z()));
    sm.add_vertex(Point_3(ma.x(), ma.y(), ma.z()));
    // put in all faces 
    vector<vector<int>> faces{ {0,6,4},{0,2,6},{0,3,2},{0,1,3},
                               {2,7,6},{2,3,7},{4,6,7},{4,7,5},
                               {0,4,5},{0,5,1},{1,5,7},{1,7,3} };
    for (const auto& f : faces) {
        std::vector<Vertex_descriptor> indices;
        for (int i : f) {
            indices.push_back(Vertex_descriptor(i));
        }
        sm.add_face(indices);
    }
    // triangualte all non-triangle faces
    PMP::triangulate_faces(sm);
    ///////////////////////////////////////////////////
    // calculate the bounding box of the compound 
    Eigen::Vector3d compBBoxMin; 
    Eigen::Vector3d compBBoxMax; 
    calculateBBox(compound, compBBoxMin, compBBoxMax);
    if (pattern.getMin().x() <= compBBoxMin.x() &&
        pattern.getMin().y() <= compBBoxMin.y() &&
        pattern.getMin().z() <= compBBoxMin.z() &&
        pattern.getMax().x() >= compBBoxMax.x() &&
        pattern.getMax().y() >= compBBoxMax.y() &&
        pattern.getMax().z() >= compBBoxMax.z()) {
        inside = compound;
        return;
    }
    ABtree tree;
    tree.add_mesh(sm, PMP::parameters::apply_per_connected_component(false));
    for (auto& com : compound.convexes) tree.add_mesh(com->convexMesh, PMP::parameters::apply_per_connected_component(false));
    std::unordered_set<int> intersected;
    for (auto& intersection : tree.get_all_intersections_and_inclusions(0)) {
        intersected.insert(intersection.first - 1);
    }
    for (int i = 0; i < compound.convexes.size(); ++i) {
        if (intersected.find(i) != intersected.end()) {
            // build inside piece first 
            auto copyMesh = compound.convexes[i]->convexMesh; 
            auto copyCell = sm;
            PMP::corefine_and_compute_intersection(copyMesh, copyCell, copyMesh);
            copyMesh.collect_garbage();
            double volume = PMP::volume(copyMesh);
            std::vector<Eigen::Vector3d> verticesInside;
            std::vector<std::vector<int>> facesInside;
            // double check that they actually intersected 
            if (volume > 0) {
                buildVFfromSM(copyMesh, verticesInside, facesInside);
                spConvex convexInside(new MeshConvex{ verticesInside, facesInside, copyMesh, {0, 0, 0}, volume });
                calculateCentroid(*convexInside, Eigen::Vector3d(0, 0, 0));
                inside.convexes.push_back(convexInside);
            }
            // build outside piece second 
            copyMesh = compound.convexes[i]->convexMesh;
            copyCell = sm;
            PMP::reverse_face_orientations(copyCell);
            PMP::corefine_and_compute_intersection(copyMesh, copyCell, copyMesh);
            copyMesh.collect_garbage();
            volume = PMP::volume(copyMesh);
            std::vector<Eigen::Vector3d> verticesOutside;
            std::vector<std::vector<int>> facesOutside;
            // double check that they actually intersected 
            if (volume > 0) {
                buildVFfromSM(copyMesh, verticesOutside, facesOutside);
                spConvex convexOutside(new MeshConvex{ verticesOutside, facesOutside, copyMesh, {0, 0, 0}, volume });
                calculateCentroid(*convexOutside, Eigen::Vector3d(0, 0, 0));
                outside.convexes.push_back(convexOutside);
            }
        }
        else {
            outside.convexes.push_back(compound.convexes[i]);
        }
    }
    inside.centroid = calculateCentroidCompound(inside.convexes);
    outside.centroid = calculateCentroidCompound(outside.convexes);
}



// The core fracture algorihm pipeline  
std::vector<Compound> fracturePipeline(Compound& compound, Pattern& pattern) {
    // pre-fracture 
    Compound inside, outside;
    preFracture(compound, pattern, inside, outside);
    // scale up each convex a little bit for intersection 
    // critical for island detection!!!!
    for (auto& c : inside.convexes) {
        for (size_t i = 0; i < c->vertices.size(); i++) {
            for (size_t i = 0; i < c->vertices.size(); i++) {
                c->vertices[i] += -c->centroid;
            }
            c->vertices[i] *= 1.01; // might need to be adjust for different mesh
            for (size_t i = 0; i < c->vertices.size(); i++) {
                c->vertices[i] += c->centroid;
            }
        }
        Surface_mesh newmesh; 
        buildSMfromVF(c->vertices, c->faces, newmesh);
        c->convexMesh = newmesh;
    }
    // scale up outside a bit as well for its island detection
    for (auto& c : outside.convexes) {
        for (size_t i = 0; i < c->vertices.size(); i++) {
            for (size_t i = 0; i < c->vertices.size(); i++) {
                c->vertices[i] += -c->centroid;
            }
            c->vertices[i] *= 1.01; // might need to be adjust for different mesh
            for (size_t i = 0; i < c->vertices.size(); i++) {
                c->vertices[i] += c->centroid;
            }
        }
        Surface_mesh newmesh;
        buildSMfromVF(c->vertices, c->faces, newmesh);
        c->convexMesh = newmesh;
    }
    //Second Step: Intersection 
    clipAABB(inside, pattern);
    //Third Step: Welding
    weldforPattern(pattern);
    //Fourth Step: Compound Formation + Island Detection
    std::vector<Compound> fractured;
    for (const auto& cell : pattern.getCells()) {
        if (cell->convexes.size() > 0) {
            auto cellCompounds = islandDetection(Compound{ cell->convexes });
            for (auto& c : cellCompounds) {
                c.centroid = calculateCentroidCompound(c.convexes);
                fractured.push_back(c);
            }
        }
    }
    auto outsideCompounds = islandDetection(outside);
    for (auto& outsidecom : outsideCompounds) {
        outsidecom.centroid = calculateCentroidCompound(outsidecom.convexes);
        fractured.push_back(outsidecom);
    }
    return fractured;
}

#if VORO_LIB
void computeVoronoiCells(
    const std::vector<vec3>& points,
    vec3 minCorner,
    vec3 maxCorner,
    AllCellVertices& cellVertices,
    AllCellFaces& cellFaces,
    AllCellEdges& cellEdges
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
#endif