/**
 * This file is to put test related functions that get called in main.cpp
**/
#pragma once
#include <iostream>

#include "clipper.h"
#include "utils.h"
#include "meshprep.h"

using namespace std;

typedef vector<vector<Eigen::Vector3d>> VoroVerts;
typedef vector<vector<vector<int>>>     VoroFaces;
typedef vector<Eigen::MatrixXi>         VoroEdges;

struct UnitCube {
    vector<Eigen::Vector3d> points{ Eigen::Vector3d(0.0,0.0,0.0),
                                            Eigen::Vector3d(0.0,0.0,1.0),
                                            Eigen::Vector3d(0.0,1.0,0.0),
                                            Eigen::Vector3d(0.0,1.0,1.0),
                                            Eigen::Vector3d(1.0,0.0,0.0),
                                            Eigen::Vector3d(1.0,0.0,1.0),
                                            Eigen::Vector3d(1.0,1.0,0.0),
                                            Eigen::Vector3d(1.0,1.0,1.0) };

    vector<vector<int>> faces{ {0,6,4},{0,2,6},{0,3,2},{0,1,3},
                                                {2,7,6},{2,3,7},{4,6,7},{4,7,5},
                                                {0,4,5},{0,5,1},{1,5,7},{1,7,3} };

    vector<Eigen::Vector3d> direction{ Eigen::Vector3d(0.25, 0.25, 0.25),
                                                Eigen::Vector3d(0.25, 0.25, -0.25),
                                                Eigen::Vector3d(0.25, -0.25, 0.25),
                                                Eigen::Vector3d(0.25, -0.25, -0.25),
                                                Eigen::Vector3d(-0.25, 0.25, 0.25),
                                                Eigen::Vector3d(-0.25, 0.25, -0.25),
                                                Eigen::Vector3d(-0.25, -0.25, 0.25),
                                                Eigen::Vector3d(-0.25, -0.25, -0.25) };
};


//////////////////////////////////////////////////////////////////////////////////////////
// TEST MESH CLIPPING
//////////////////////////////////////////////////////////////////////////////////////////

MeshConvex testFunc2(int cellIndex, 
    const Eigen::MatrixXd& V, const Eigen::MatrixXi& F,
    const VoroVerts& cellVertices,
    const VoroFaces& cellFaces,
    const VoroEdges& cellEdges
) {
    /*std::vector<Eigen::Vector3d> points = { Eigen::Vector3d(0.0,0.0,0.0),
                                            Eigen::Vector3d(0.0,0.0,1.0),
                                            Eigen::Vector3d(0.0,1.0,0.0),
                                            Eigen::Vector3d(0.0,1.0,1.0),
                                            Eigen::Vector3d(1.0,0.0,0.0),
                                            Eigen::Vector3d(1.0,0.0,1.0),
                                            Eigen::Vector3d(1.0,1.0,0.0),
                                            Eigen::Vector3d(1.0,1.0,1.0) };
    for (auto& p : points) {
        p += Eigen::Vector3d(-0.5, -0.5, -0.5);
        p *= 1;
    }
    std::vector<std::vector<int>> faces = { {0,6,4},{0,2,6},{0,3,2},{0,1,3},
                                            {2,7,6},{2,3,7},{4,6,7},{4,7,5},
                                            {0,4,5},{0,5,1},{1,5,7},{1,7,3} };*/
    
    // hard-coded geometry pre-processing
    std::vector<Eigen::Vector3d> points;
    convertToVertArray(V, points);
    
    std::vector<std::vector<int>> faces;
    convertToFaceVertArray(F, faces);
    
    // convert geometry to a CGAL surface mesh
    Surface_mesh tester_mesh;
    buildSMfromVF(points, faces, tester_mesh);
    
    // generate 3d voronoi diagram
    Pattern pattern(cellVertices, cellFaces, cellEdges);
    pattern.createCellsfromVoro();
    Cell cell = *pattern.getCells()[cellIndex];
    
    MeshConvex tester{ points, faces, tester_mesh };
    spConvex result(new MeshConvex);
    clipConvexAgainstCell(tester, cell, result);
    
    // clipped convex post-processing for visualization
    calculateCentroid(*result, Eigen::Vector3d(0, 0, 0));
    translateMesh(*result, result->centroid, 0.0);
    
    return *result;
}

// Compute every clipped convex in a voronoi cell 
MeshConvex testFunc3(const Eigen::MatrixXd& V, const Eigen::MatrixXi& F,
    const VoroVerts& cellVertices,
    const VoroFaces& cellFaces,
    const VoroEdges& cellEdges
) {
    // Convert format from V,F matrices
    std::vector<Eigen::Vector3d> points;
    convertToVertArray(V, points);
    std::vector<std::vector<int>> faces;
    convertToFaceVertArray(F, faces);
    
    // Build SM of the whole mesh
    Surface_mesh tester_mesh;
    buildSMfromVF(points, faces, tester_mesh);
    
    // Calculate the centroid of the whole mesh use for visualization
    MeshConvex tester{ points, faces, tester_mesh };
    calculateCentroid(tester, Eigen::Vector3d(0, 0, 0));
    
    // Construct Pattern from Voronoi Decomposition
    Pattern pattern(cellVertices, cellFaces, cellEdges);
    pattern.createCellsfromVoro();
    auto cells = pattern.getCells();
    
    std::vector<Eigen::Vector3d> final_vertices;
    std::vector<std::vector<int>> final_faces;
    
    // Calculate the center of mass for each clipped mesh and 
    // move it away from the center of mass of the whole mesh
    for (auto const& c : cells) {
        int previous_verts = final_vertices.size();
        spConvex result(new MeshConvex);
        clipConvexAgainstCell(tester, *c, result);
        
        calculateCentroid(*result, tester.centroid);
        translateMesh(*result, result->centroid, 0.03);
        
        for (size_t i = 0; i < result->faces.size(); i++) {
            result->faces[i][0] += previous_verts;
            result->faces[i][1] += previous_verts;
            result->faces[i][2] += previous_verts;
        }
        final_vertices.insert(final_vertices.end(), result->vertices.begin(), result->vertices.end());
        final_faces.insert(final_faces.end(), result->faces.begin(), result->faces.end());
    }

    return MeshConvex{ final_vertices, final_faces };
}

//////////////////////////////////////////////////////////////////////////////////////////
// TEST WELDING
//////////////////////////////////////////////////////////////////////////////////////////

void pre_test_welding(Eigen::MatrixXd& V, Eigen::MatrixXi& F) {
    UnitCube cube;
    for (auto& p : cube.points) {
        p += Eigen::Vector3d(-0.5, -0.5, -0.5);
        //p *= 1;
    }

    V = convertToMatrixXd(cube.points);
    F = convertToMatrixXi(cube.faces);
}

//Functional calls to test welding mechanism, it manually creates 8 smaller 
//individual cubes forming a single large cube.
void testWelding(std::vector<MeshConvex>& clippedMeshConvex,
    const VoroVerts& cellVertices,
    const VoroFaces& cellFaces,
    const VoroEdges& cellEdges,
    std::vector<Pattern::sPCell>& gCells
) {
    UnitCube cube;

    // move each small cube to their correct position and scaling
    std::vector<MeshConvex> allCubes;
    for (size_t i = 0; i < 8; ++i) {
        auto p = cube.points;
        for (auto& eachp : p) {
            eachp += Eigen::Vector3d(-0.5, -0.5, -0.5);
            eachp *= 0.5;
            eachp += cube.direction[i];
        }
        Surface_mesh sm;
        buildSMfromVF(p, cube.faces, sm);
        allCubes.push_back(MeshConvex{ p, cube.faces, sm });
    }
    
    clippedMeshConvex = allCubes; // for testing the cube
    
    // generate voronoi diagram
    Pattern pattern(cellVertices, cellFaces, cellEdges);
    pattern.createCellsfromVoro();
    
    for (auto& c : pattern.getCells()) {
        for (auto& cube : allCubes) {
            spConvex clipped(new MeshConvex);
            clipConvexAgainstCell(cube, *c, clipped);
            // only add to cell's list if intersected
            if (clipped->volume > 0) {
                c->convexes.push_back(clipped);
            }
        }
    }
    
    // welding
    weldforPattern(pattern);
    double sum_v = 0;
    gCells = pattern.getCells();
    for (const auto& c : gCells) {
        sum_v += c->convexes[0]->volume;
    }
    std::cout << "Total volume of summing each Cell's first convex piece: " << sum_v << std::endl;
}

//////////////////////////////////////////////////////////////////////////////////////////
// TESTING HASHING Eigen::Vector4d
//////////////////////////////////////////////////////////////////////////////////////////

void testHashing() {
    std::vector<Eigen::Vector3d> t_points;
    // generate 3 random points
    for (size_t i = 0; i < 3; i++) {
        t_points.push_back(Eigen::Vector3d(double(rand()) / RAND_MAX, double(rand()) / RAND_MAX, double(rand()) / RAND_MAX));
    }
    // creating 2 opposite planes 
    Eigen::Vector4d plane1 = createPlane(t_points[0], t_points[1], t_points[2]);
    Eigen::Vector4d plane2 = createPlane(t_points[2], t_points[1], t_points[0]);
    std::cout << "The first plane:\n" << plane1 << std::endl;
    std::cout << "The second plane:\n" << plane2 << std::endl;
    auto plane_map = std::unordered_map<Eigen::Vector4d, std::string>();
    plane_map[plane1] = "hhh";
    plane_map[plane2] = "bbbb";
    bool checker = plane_map.find(plane1) != plane_map.end();
    bool checker2 = plane_map.find(plane2) != plane_map.end();
    bool checker3 = plane_map.find(Eigen::Vector4d(0.1, 0.2, 0.3, 0.4)) != plane_map.end();
}

//////////////////////////////////////////////////////////////////////////////////////
// TESTING ISLAND DETECTION
//////////////////////////////////////////////////////////////////////////////////////////

void testIsland(std::vector<MeshConvex>& clippedMeshConvex,
    const VoroVerts& cellVertices,
    const VoroFaces& cellFaces,
    const VoroEdges& cellEdges,
    std::vector<Compound>& gCompounds,
    bool isCustomMesh = false
) { 
    if (!isCustomMesh) {
        UnitCube cube;
        std::vector<Eigen::Vector3d> direction{ Eigen::Vector3d(0.25, 0.25, 0.25),
                                                Eigen::Vector3d(0.25, 0.25, -0.25),
            //Eigen::Vector3d(0.25, -0.25, 0.25),
            //Eigen::Vector3d(0.25, -0.25, -0.25),
            //Eigen::Vector3d(-0.25, 0.25, 0.25),
            //Eigen::Vector3d(-0.25, 0.25, -0.25),
            Eigen::Vector3d(-0.25, -0.25, 0.25),
            Eigen::Vector3d(-0.25, -0.25, -0.25)
        };

        // move each small cube to their correct position and scaling
        std::vector<MeshConvex> allCubes;
        for (size_t i = 0; i < direction.size(); ++i) {
            auto p = cube.points;
            for (auto& eachp : p) {
                eachp += Eigen::Vector3d(-0.5, -0.5, -0.5);
                eachp *= 0.5;
                eachp += direction[i];
            }
            Surface_mesh sm;
            buildSMfromVF(p, cube.faces, sm);
            allCubes.push_back(MeshConvex{ p, cube.faces, sm });
        }

        clippedMeshConvex = allCubes; // for testing the cube
    }
    
    // generate voronoi diagram
    Pattern pattern(cellVertices, cellFaces, cellEdges);
    pattern.createCellsfromVoro();

    for (auto& cell : pattern.getCells()) {
        for (auto& convex : clippedMeshConvex) {
            spConvex clipped(new MeshConvex);
            clipConvexAgainstCell(convex, *cell, clipped);
            // only add to cell's list if intersected
            if (clipped->volume > 0) {
                cell->convexes.push_back(clipped);
            }
        }
    }
    // compound formation 
    for (const auto& cell : pattern.getCells()) {
        auto com = Compound{ cell->convexes };
        gCompounds.push_back(com);
    }
}

////////////////////////////////////////////////////////////////////
//TEST for Customize Readobj function
////////////////////////////////////////////////////////////////////
void testObj(const std::string& filePath, 
    std::vector<spConvex>& convexes) {
    convexes = mesh_prep::readOBJ(filePath);
}