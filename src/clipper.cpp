#include "clipper.h"
#include <iostream>
#include <chrono>

#define MULTITHREAD 
#ifdef MULTITHREAD
    #define CGAL_HAS_THREADS
    #define BOOST_HAS_THREADS
#endif
std::mutex patternMutex;
Pattern::Pattern(AllCellVertices v, AllCellFaces f, AllCellEdges e):o_cellVertices(v), o_cellFaces(f), o_cellEdges(e)
{}

std::vector<Pattern::spCell> Pattern::getCells() {
    return cells;
}

Pattern::AllCellVertices Pattern::getVertices() {
    return o_cellVertices;
}

Pattern::AllCellFaces Pattern::getFaces() {
    return o_cellFaces;
}

void Pattern::createCellsfromVoro() {
    // for every cell create Cell and push to container
    for (int i = 0; i < o_cellFaces.size(); ++i) {
        auto curCell = o_cellFaces[i];
        Surface_mesh sm; 
        buildSMfromVF(o_cellVertices[i], o_cellFaces[i], sm);
        spCell cellptr(new Cell{i, std::vector<spConvex>{}, o_cellVertices[i], o_cellFaces[i], sm });
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

void translateMesh(MeshConvex& mesh, Eigen::Vector3d direction, double scale) {
    auto d = direction.normalized();
    for (size_t i = 0; i < mesh.vertices.size(); i++) {
        mesh.vertices[i] += d * scale;
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

#ifdef MULTITHREAD
void workerClipAABB(int cellID, int cell_num, Compound& compound, Pattern& pattern, std::vector<std::pair<size_t, bool>> intersects) {
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
            }
            spConvex result(new MeshConvex{ vertices, faces, convexm, {0, 0, 0}, volume });
            pattern.getCells()[cellID]->convexes.push_back(result);
        }
    }
}
#endif

void clipAABB(Compound& compound, Pattern& pattern) {
    ABtree tree;
    int cell_num = pattern.getCells().size();
    int convex_num = compound.convexes.size();
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
#ifdef MULTITHREAD
    std::vector<std::thread> threads;
    for (int i = 0; i < cell_num; ++i) {
        std::thread myThread(workerClipAABB, i, cell_num, compound, pattern, intersection_list[i]);
        threads.push_back(std::move(myThread));
    }
    for (auto& thread : threads) {
        thread.join();
    }
#endif
#ifndef MULTITHREAD
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
        if (total_volume - sum_volume < DBL_EPSILON) {
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

// A custom intersection function used in island detection 
template<typename T>
std::unordered_set<T> intersection(const std::unordered_set<T>& set1, const std::unordered_set<T>& set2) {
    std::unordered_set<T> result;

    // Iterate through the first set and check if each element is in the second set
    for (const T& element : set1) {
        if (set2.find(element) != set2.end()) {
            // If the element is found in the second set, add it to the result
            result.insert(element);
        }
    }

    return result;
}

std::vector<Compound> islandDetection(Compound& old_compound) {
    // store mapping relation between a face plane and the group index
    std::unordered_map<Eigen::Vector4d, std::shared_ptr<MeshConvex>> linker;
    std::vector<Compound> results; 
    std::vector<std::shared_ptr<MeshConvex>> convs = old_compound.convexes;
    // No need for island detection if there's only 1 piece
    if (convs.size() < 2) return std::vector<Compound>{old_compound};
    // clear previous groupings 
    for (auto& c : old_compound.convexes) {
        c->group.clear();
    }
    // iterate over every face of every convex piece
    for (int i = 0; i < convs.size(); ++i) {
        auto convex = convs[i];
        auto vertices = convex->vertices;
        auto cur = std::unordered_set<int>{ i };
        std::vector<Eigen::Vector4d> planes; 
        // first iteration update the set of convexes
        for (const auto& face : convex->faces) {
            auto p1 = vertices[face[0]];
            auto p2 = vertices[face[1]];
            auto p3 = vertices[face[2]];
            auto plane = createPlane(p1, p2, p3);
            planes.push_back(plane);
            if (linker.find(-plane) != linker.end()) {
                auto prev = linker[-plane]->group;
                cur.insert(prev.begin(), prev.end());
            }
        }
        // second iteration update to group variable of convexes
        for (const auto& p : planes) {
            linker[p] = convs[i];
            linker[p]->group = cur;
            if (linker.find(-p) != linker.end()) {
                std::vector<int> temp_group; 
                for (const auto& prev_group : linker[-p]->group) convs[prev_group]->group = cur;
                linker[-p]->group = cur;
            }
        }
    }
    // use convs to group islands together
    std::unordered_set<int> cur_convs;
    std::unordered_set<int> final_convs;
    for (int l = 0; l < convs.size(); ++l) final_convs.insert(l);
    for (const auto& c : convs) {
        if (cur_convs == final_convs) break;
        auto g = c->group;
        std::unordered_set<int> inter = intersection(g, cur_convs);
        //std::set_intersection(cur_convs.begin(), cur_convs.end(), g.begin(), g.end(),
        //    std::inserter(inter, inter.begin()));
        if (inter.size() == 0) {
            cur_convs.insert(g.begin(), g.end());
            std::vector<std::shared_ptr<MeshConvex>> new_convex;
            for (const auto& index : g) new_convex.push_back(convs[index]);
            results.push_back(Compound{ new_convex });
        }
    }
    return results;
}
// alternative island detection using intersection 
// not correct implementation
std::vector<Compound> islandDetection2(Compound& old_compound) {
    std::unordered_set<int> cur_convs;
    std::unordered_set<int> final_convs;
    std::vector<Compound> results; 
    int N = old_compound.convexes.size();
    for (int l = 0; l < N; ++l) final_convs.insert(l);
    for (int i = 0; i < N; ++i) {
        if (cur_convs == final_convs) break;
        std::unordered_set<int> curr; 
        if (cur_convs.find(i) == cur_convs.end()) {
            curr.insert(i);
            for (int j = i+1; j < N; ++j) {
                if (PMP::do_intersect(old_compound.convexes[i]->convexMesh, old_compound.convexes[j]->convexMesh,
                    PMP::parameters::do_overlap_test_of_bounded_sides(true), PMP::parameters::do_overlap_test_of_bounded_sides(true))) {
                    curr.insert(j);
                }
            }
            std::vector<spConvex> new_convex;
            for (const auto& index : curr) {
                new_convex.push_back(old_compound.convexes[index]);
                cur_convs.insert(index);
            }
            results.push_back(Compound{ new_convex });
        }
    }
    return results;
}


// The core fracture algorihm pipeline  
std::vector<Compound> fracturePipeline(Compound& compound, Pattern& pattern) {
    //TODO: alignmnet First Step: Alignment
    //TODO: acceleration structure!!!!! Second Step: Intersection 
    /*for (auto& cell : pattern.getCells()) {
        for (auto& convex : compound.convexes) {
            spConvex clipped(new MeshConvex);
            if (clipConvexAgainstCell(*convex, *cell, clipped)) cell->convexes.push_back(clipped);
        }
    }*/
    clipAABB(compound, pattern);
    //Third Step: Welding
    weldforPattern(pattern);
    //Fourth Step: Compound Formation + Island Detection
    std::vector<Compound> fractured;
    for (const auto& cell : pattern.getCells()) {
        if (cell->convexes.size() > 0) {
            // DISABLE island detection for now 
            /*auto cellCompounds = islandDetection2(Compound{ cell->convexes });
            for (const auto& c : cellCompounds) fractured.push_back(c);*/
            fractured.push_back(Compound{ cell->convexes });
        }
    }
    return fractured;
}
