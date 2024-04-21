/**
 * This file is to put test related functions that get called in main.cpp
**/
#pragma once
#include <iostream>

#include "utils.h"
#include "meshprep.h"

using namespace std;

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

// Input polygon
Eigen::MatrixXd V;                              // #V by 3 matrix for vertices
Eigen::MatrixXi F;                              // matrix for face indices
Eigen::MatrixXd B;                              // matrix for barycenters
Eigen::MatrixXd N;                              // matrix for normals
Eigen::Vector3d minCorner, maxCorner;           // min and max corners of mesh's bounding box

Eigen::MatrixXd meshV;
Eigen::MatrixXi meshF;

// Tetrahedralized interior
Eigen::MatrixXd TV;                             // #TV by 3 matrix for vertex positions
Eigen::MatrixXi TT;                             // #TT by 4 matrix for tet face indices
Eigen::MatrixXi TF;                             // #TF by 3 matrix for triangle face indices ('f', else `boundary_facets` is called on TT)

// Voronoi diagram
int gNumPoints = 10;
vector<Eigen::Vector3d> gPoints =               // cell nodes
{
  Eigen::Vector3d(-0.0695318, -0.0872571, 0.215309),
  Eigen::Vector3d(-0.242122, 0.175633, -0.459787),
  Eigen::Vector3d(-0.382464, -0.491852, 0.480495),
  Eigen::Vector3d(0.103244, 0.349015, 0.0440263),
  Eigen::Vector3d(-0.448226, 0.216952, -0.236705),
  Eigen::Vector3d(0.097456, 0.367117, -0.0303783),
  Eigen::Vector3d(-0.318482, 0.0231078, 0.142635),
  Eigen::Vector3d(0.0593644, 0.111092, 0.203777),
  Eigen::Vector3d(-0.356605, 0.474539, -0.167495),
  Eigen::Vector3d(0.410755, 0.111101, 0.279034)
};
vector<vec3> gPointsVec;

// hard-coded 10 Voronoi cells for now
AllCellVertices gCellVertices =   // Vertices of each Voronoi cell
{
  {
      Eigen::Vector3d(-0.5, -0.5, -0.5),
      Eigen::Vector3d(-0.0740757, 0.0499249, -0.144954),
      Eigen::Vector3d(-0.5, 0.243325, -0.5),
      Eigen::Vector3d(-0.0245021, 0.5, -0.384857),
      Eigen::Vector3d(-0.128199, 0.12096, -0.13383),
      Eigen::Vector3d(-0.0695257, 0.032968, -0.148671),
      Eigen::Vector3d(-0.318479, 0.5, -0.5),
      Eigen::Vector3d(-0.147488, 0.247829, -0.17515),
      Eigen::Vector3d(-0.5, -0.5, -0.36232),
      Eigen::Vector3d(-0.469712, -0.5, -0.334337),
      Eigen::Vector3d(-0.459359, 0.446042, -0.5),
      Eigen::Vector3d(-0.259384, -0.5, -0.307677),
      Eigen::Vector3d(0.5, -0.171943, -0.5),
      Eigen::Vector3d(0.1211, 0.5, -0.5),
      Eigen::Vector3d(0.5, -0.199609, -0.487663),
      Eigen::Vector3d(0.366353, -0.222884, -0.371595),
      Eigen::Vector3d(0.367306, -0.264124, -0.376039),
      Eigen::Vector3d(0.480758, -0.5, -0.496896),
      Eigen::Vector3d(0.48427, -0.5, -0.5),
      Eigen::Vector3d(0.5, -0.340857, -0.5)
  },
  {
      Eigen::Vector3d(-0.459359, 0.446042, -0.5),
      Eigen::Vector3d(-0.469712, -0.5, -0.334337),
      Eigen::Vector3d(-0.5, 0.460498, -0.5),
      Eigen::Vector3d(-0.5, 0.243325, -0.5),
      Eigen::Vector3d(-0.5, -0.5, -0.323977),
      Eigen::Vector3d(-0.5, -0.5, -0.36232),
      Eigen::Vector3d(-0.5, 0.302956, 0.0863368),
      Eigen::Vector3d(-0.128199, 0.12096, -0.13383),
      Eigen::Vector3d(-0.147488, 0.247829, -0.17515),
      Eigen::Vector3d(-0.187059, 0.23037, -0.0577889)
  },
  {
      Eigen::Vector3d(-0.5, 0.460498, -0.5),
      Eigen::Vector3d(-0.0245021, 0.5, -0.384857),
      Eigen::Vector3d(-0.5, 0.5, -0.5),
      Eigen::Vector3d(-0.147488, 0.247829, -0.17515),
      Eigen::Vector3d(-0.5, 0.302956, 0.0863368),
      Eigen::Vector3d(-0.187059, 0.23037, -0.0577889),
      Eigen::Vector3d(-0.5, 0.5, 0.373158),
      Eigen::Vector3d(-0.159309, 0.5, 0.0615546),
      Eigen::Vector3d(-0.459359, 0.446042, -0.5),
      Eigen::Vector3d(-0.318479, 0.5, -0.5),
      Eigen::Vector3d(-0.290812, 0.5, 0.347443),
      Eigen::Vector3d(-0.196704, 0.276124, 0.0099971)
  },
  {
      Eigen::Vector3d(0.366353, -0.222884, -0.371595),
      Eigen::Vector3d(0.5, 0.361425, -0.0234486),
      Eigen::Vector3d(0.1211, 0.5, -0.5),
      Eigen::Vector3d(0.5, 0.5, -0.5),
      Eigen::Vector3d(-0.159309, 0.5, 0.0615546),
      Eigen::Vector3d(0.304808, 0.0490665, -0.0842575),
      Eigen::Vector3d(-0.0245021, 0.5, -0.384857),
      Eigen::Vector3d(0.5, 0.5, 0.0102653),
      Eigen::Vector3d(-0.12193, 0.169778, -0.0216926),
      Eigen::Vector3d(-0.196704, 0.276124, 0.0099971),
      Eigen::Vector3d(-0.0740757, 0.0499249, -0.144954),
      Eigen::Vector3d(0.5, -0.199609, -0.487663),
      Eigen::Vector3d(-0.187059, 0.23037, -0.0577889),
      Eigen::Vector3d(0.5, -0.171943, -0.5),
      Eigen::Vector3d(-0.147488, 0.247829, -0.17515),
      Eigen::Vector3d(-0.128199, 0.12096, -0.13383)
  },
  {
      Eigen::Vector3d(-0.5, -0.5, -0.323977),
      Eigen::Vector3d(-0.0695257, 0.032968, -0.148671),
      Eigen::Vector3d(-0.271749, 0.5, 0.428967),
      Eigen::Vector3d(-0.283244, 0.5, 0.5),
      Eigen::Vector3d(-0.324052, -0.114024, 0.5),
      Eigen::Vector3d(-0.5, -0.0921634, 0.5),
      Eigen::Vector3d(-0.5, 0.5, 0.5),
      Eigen::Vector3d(-0.203555, 0.157781, 0.5),
      Eigen::Vector3d(-0.259384, -0.5, -0.307677),
      Eigen::Vector3d(-0.5, 0.5, 0.373158),
      Eigen::Vector3d(-0.5, -0.5, -0.121618),
      Eigen::Vector3d(-0.323459, -0.5, -0.0881849),
      Eigen::Vector3d(-0.290812, 0.5, 0.347443),
      Eigen::Vector3d(-0.5, 0.302956, 0.0863368),
      Eigen::Vector3d(-0.469712, -0.5, -0.334337),
      Eigen::Vector3d(-0.12193, 0.169778, -0.0216926),
      Eigen::Vector3d(-0.187059, 0.23037, -0.0577889),
      Eigen::Vector3d(-0.196704, 0.276124, 0.0099971),
      Eigen::Vector3d(-0.0740757, 0.0499249, -0.144954),
      Eigen::Vector3d(-0.128199, 0.12096, -0.13383)
  },
  {
      Eigen::Vector3d(0.304808, 0.0490665, -0.0842575),
      Eigen::Vector3d(0.5, 0.5, 0.11684),
      Eigen::Vector3d(-0.271749, 0.5, 0.428967),
      Eigen::Vector3d(0.5, 0.5, 0.0102653),
      Eigen::Vector3d(-0.0131432, 0.5, 0.5),
      Eigen::Vector3d(0.20718, 0.5, 0.5),
      Eigen::Vector3d(-0.12193, 0.169778, -0.0216926),
      Eigen::Vector3d(0.179669, 0.46444, 0.5),
      Eigen::Vector3d(0.5, 0.361425, -0.0234486),
      Eigen::Vector3d(-0.159309, 0.5, 0.0615546),
      Eigen::Vector3d(-0.196704, 0.276124, 0.0099971),
      Eigen::Vector3d(-0.290812, 0.5, 0.347443)
  },
  {
      Eigen::Vector3d(0.480758, -0.5, -0.496896),
      Eigen::Vector3d(0.367306, -0.264124, -0.376039),
      Eigen::Vector3d(0.174981, -0.5, 0.5),
      Eigen::Vector3d(-0.0695257, 0.032968, -0.148671),
      Eigen::Vector3d(-0.324052, -0.114024, 0.5),
      Eigen::Vector3d(0.34849, -0.5, 0.5),
      Eigen::Vector3d(-0.203555, 0.157781, 0.5),
      Eigen::Vector3d(0.179682, -0.0912641, 0.5),
      Eigen::Vector3d(-0.259384, -0.5, -0.307677),
      Eigen::Vector3d(-0.323459, -0.5, -0.0881849)
  },
  {
      Eigen::Vector3d(-0.271749, 0.5, 0.428967),
      Eigen::Vector3d(0.367306, -0.264124, -0.376039),
      Eigen::Vector3d(-0.0131432, 0.5, 0.5),
      Eigen::Vector3d(0.179669, 0.46444, 0.5),
      Eigen::Vector3d(-0.203555, 0.157781, 0.5),
      Eigen::Vector3d(0.179682, -0.0912641, 0.5),
      Eigen::Vector3d(-0.283244, 0.5, 0.5),
      Eigen::Vector3d(-0.0695257, 0.032968, -0.148671),
      Eigen::Vector3d(0.304808, 0.0490665, -0.0842575),
      Eigen::Vector3d(-0.12193, 0.169778, -0.0216926),
      Eigen::Vector3d(-0.0740757, 0.0499249, -0.144954),
      Eigen::Vector3d(0.366353, -0.222884, -0.371595)
  },
  {
      Eigen::Vector3d(0.34849, -0.5, 0.5),
      Eigen::Vector3d(0.5, -0.5, -0.5),
      Eigen::Vector3d(0.20718, 0.5, 0.5),
      Eigen::Vector3d(0.48427, -0.5, -0.5),
      Eigen::Vector3d(0.179682, -0.0912641, 0.5),
      Eigen::Vector3d(0.5, -0.5, 0.5),
      Eigen::Vector3d(0.5, 0.5, 0.11684),
      Eigen::Vector3d(0.5, 0.5, 0.5),
      Eigen::Vector3d(0.304808, 0.0490665, -0.0842575),
      Eigen::Vector3d(0.366353, -0.222884, -0.371595),
      Eigen::Vector3d(0.5, 0.361425, -0.0234486),
      Eigen::Vector3d(0.179669, 0.46444, 0.5),
      Eigen::Vector3d(0.5, -0.199609, -0.487663),
      Eigen::Vector3d(0.5, -0.340857, -0.5),
      Eigen::Vector3d(0.480758, -0.5, -0.496896),
      Eigen::Vector3d(0.367306, -0.264124, -0.376039)
  },
  {
      Eigen::Vector3d(-0.5, -0.5, 0.5),
      Eigen::Vector3d(-0.323459, -0.5, -0.0881849),
      Eigen::Vector3d(-0.5, -0.5, -0.121618),
      Eigen::Vector3d(0.174981, -0.5, 0.5),
      Eigen::Vector3d(-0.5, -0.0921634, 0.5),
      Eigen::Vector3d(-0.324052, -0.114024, 0.5)
  }
};
AllCellFaces gCellFaces =           // Faces of each Voronoi cell represented by vertex indices
{
  {
      {5, 16, 15, 1},
      {15, 14, 12, 13, 3, 7, 4, 1},
      {4, 9, 11, 5, 1},
      {8, 9, 4, 7, 10, 2},
      {10, 6, 13, 12, 19, 18, 0, 2},
      {0, 8, 2},
      {13, 6, 3},
      {6, 10, 7, 3},
      {11, 17, 16, 5},
      {0, 18, 17, 11, 9, 8},
      {14, 19, 12},
      {15, 16, 17, 18, 19, 14}
  },
  {
      {7, 9, 6, 4, 1},
      {4, 5, 1},
      {5, 3, 0, 8, 7, 1},
      {6, 9, 8, 0, 2},
      {0, 3, 2},
      {3, 5, 4, 6, 2},
      {8, 9, 7}
  },
  {
      {7, 11, 5, 3, 1},
      {3, 8, 9, 1},
      {9, 2, 6, 10, 7, 1},
      {9, 8, 0, 2},
      {0, 4, 6, 2},
      {5, 4, 0, 8, 3},
      {5, 11, 10, 6, 4},
      {10, 11, 7}
  },
  {
      {5, 0, 11, 1},
      {11, 13, 3, 7, 1},
      {7, 4, 9, 8, 5, 1},
      {13, 11, 0, 10, 15, 14, 6, 2},
      {6, 4, 7, 3, 2},
      {3, 13, 2},
      {6, 14, 12, 9, 4},
      {8, 10, 0, 5},
      {9, 12, 15, 10, 8},
      {14, 15, 12}
  },
  {
      {8, 14, 19, 18, 1},
      {18, 15, 2, 3, 7, 1},
      {7, 4, 11, 8, 1},
      {15, 17, 12, 2},
      {12, 9, 6, 3, 2},
      {6, 5, 4, 7, 3},
      {5, 10, 11, 4},
      {6, 9, 13, 0, 10, 5},
      {11, 10, 0, 14, 8},
      {12, 17, 16, 13, 9},
      {16, 19, 14, 0, 13},
      {18, 19, 16, 17, 15}
  },
  {
      {5, 7, 0, 8, 1},
      {8, 3, 1},
      {3, 9, 11, 2, 4, 5, 1},
      {11, 10, 6, 2},
      {6, 0, 7, 4, 2},
      {8, 0, 6, 10, 9, 3},
      {7, 5, 4},
      {10, 11, 9}
  },
  {
      {0, 8, 3, 1},
      {3, 6, 7, 1},
      {7, 5, 0, 1},
      {4, 9, 2},
      {9, 8, 0, 5, 2},
      {5, 7, 6, 4, 2},
      {8, 9, 4, 6, 3}
  },
  {
      {7, 10, 11, 1},
      {11, 8, 3, 5, 1},
      {5, 4, 7, 1},
      {3, 8, 9, 0, 2},
      {0, 6, 2},
      {6, 4, 5, 3, 2},
      {6, 0, 9, 10, 7, 4},
      {11, 10, 9, 8}
  },
  {
      {3, 13, 1},
      {13, 12, 10, 6, 7, 5, 1},
      {5, 0, 14, 3, 1},
      {6, 10, 8, 11, 2},
      {11, 4, 0, 5, 7, 2},
      {7, 6, 2},
      {14, 15, 9, 12, 13, 3},
      {15, 14, 0, 4},
      {11, 8, 9, 15, 4},
      {10, 12, 9, 8}
  },
  {
      {2, 4, 5, 1},
      {5, 3, 1},
      {3, 0, 2, 1},
      {0, 4, 2},
      {5, 4, 0, 3}
  }
};
// TODO: vert indices order below is incorrect, but since we're not using this info, just leave it here to avoid visualization error. 
AllCellEdges gCellEdges =              // Edges of each Voronoi cell represented by vertex indices
{
  (Eigen::MatrixXi(30, 2) << 1, 15, 5, 1, 1, 4, 3, 3, 12, 12, 14, 5, 9, 4, 2, 7, 8, 2, 0, 0, 18, 12, 6, 6, 0, 3, 16, 11, 17, 14, 15, 16, 16, 5, 4, 7, 7, 13, 13, 14, 15, 11, 11, 9, 10, 10, 9, 8, 2, 18, 19, 19, 13, 10, 8, 6, 17, 17, 18, 19).finished(),
  (Eigen::MatrixXi(15, 2) << 1, 4, 6, 7, 1, 1, 4, 7, 0, 0, 3, 0, 8, 2, 2, 4, 6, 9, 9, 7, 5, 5, 8, 8, 3, 5, 2, 9, 6, 3).finished(),
  (Eigen::MatrixXi(18, 2) << 1, 3, 5, 7, 1, 1, 8, 3, 7, 6, 2, 2, 0, 0, 4, 0, 4, 10, 3, 5, 11, 11, 7, 9, 9, 8, 10, 10, 6, 9, 2, 8, 6, 4, 5, 11).finished(),
  (Eigen::MatrixXi(24, 2) << 1, 0, 0, 1, 1, 3, 3, 11, 5, 8, 4, 4, 2, 6, 14, 10, 0, 2, 2, 4, 9, 12, 8, 12, 11, 11, 5, 5, 7, 7, 13, 13, 8, 9, 9, 7, 6, 14, 15, 15, 10, 13, 3, 6, 12, 14, 10, 15).finished(),
  (Eigen::MatrixXi(30, 2) << 1, 18, 14, 8, 1, 1, 3, 2, 2, 15, 8, 4, 4, 2, 12, 15, 3, 6, 9, 4, 5, 10, 5, 0, 0, 9, 0, 13, 16, 16, 18, 19, 19, 14, 8, 7, 7, 3, 15, 18, 11, 11, 7, 12, 17, 17, 6, 9, 12, 5, 6, 11, 10, 10, 13, 13, 14, 16, 17, 19).finished(),
  (Eigen::MatrixXi(18, 2) << 1, 0, 0, 5, 1, 1, 3, 4, 2, 2, 9, 3, 2, 6, 10, 4, 0, 9, 8, 8, 7, 7, 5, 3, 8, 5, 4, 11, 11, 9, 6, 10, 11, 7, 6, 10).finished(),
  (Eigen::MatrixXi(15, 2) << 1, 3, 0, 0, 1, 6, 3, 0, 5, 2, 4, 2, 2, 8, 4, 3, 8, 8, 1, 7, 7, 6, 5, 7, 9, 9, 4, 5, 9, 6).finished(),
  (Eigen::MatrixXi(18, 2) << 1, 10, 7, 1, 1, 3, 3, 8, 4, 4, 0, 0, 8, 2, 2, 0, 4, 9, 11, 11, 10, 7, 5, 5, 8, 11, 7, 5, 2, 9, 9, 3, 6, 6, 6, 10).finished(),
  (Eigen::MatrixXi(24, 2) << 1, 3, 1, 1, 5, 6, 6, 10, 12, 3, 0, 0, 2, 8, 8, 2, 2, 0, 4, 9, 9, 14, 4, 8, 13, 13, 3, 5, 7, 7, 10, 12, 13, 14, 14, 5, 11, 11, 10, 6, 7, 4, 11, 12, 15, 15, 15, 9).finished(),
  (Eigen::MatrixXi(9, 2) << 1, 4, 2, 1, 1, 3, 0, 0, 0, 5, 5, 4, 2, 3, 5, 2, 3, 4).finished()
};

// Mesh operations
std::vector<MeshConvex> gClippedMeshConvex;     // Global var for testing mesh clipping 
std::vector<spCell> gCells;            // Global var for testing welding 
std::vector<Compound> gCompounds;               // Global var for testing island detection 
std::vector<Compound> gCurrCompounds;           // Global var for testing island detection 
int  gCurrConvex = 0;                           // Global var for testing island detection
std::vector<spConvex> gFracturedConvex;         // Global var for testing pipeline

// ReadObj testing 
Compound ginitialConvexes;                      // Global var for testing readOBJ function
std::string gOBJPath = "..\\assets\\results\\bunny_out.obj";


void generateRandomPoints(int numPoints, std::vector<Eigen::Vector3d>& points)
{
    points.clear();

    std::random_device rd;  // Obtain random number from hardware and seed the generator
    std::mt19937 gen(rd());
    //std::mt19937 gen(19);
    std::uniform_real_distribution<> disX(minCorner.x(), maxCorner.x());
    std::uniform_real_distribution<> disY(minCorner.y(), maxCorner.y());
    std::uniform_real_distribution<> disZ(minCorner.z(), maxCorner.z());

    for (int i = 0; i < numPoints; ++i) {
        double x = disX(gen);
        double y = disY(gen);
        double z = disZ(gen);
        points.push_back(Eigen::Vector3d(x, y, z));
    }
}

//////////////////////////////////////////////////////////////////////////////////////////
// TEST MESH CLIPPING
//////////////////////////////////////////////////////////////////////////////////////////

MeshConvex testFunc2(int cellIndex, 
    const Eigen::MatrixXd& V, const Eigen::MatrixXi& F,
    const AllCellVertices& cellVertices,
    const AllCellFaces& cellFaces,
    const AllCellEdges& cellEdges
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
    const AllCellVertices& cellVertices,
    const AllCellFaces& cellFaces,
    const AllCellEdges& cellEdges
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
    const AllCellVertices& cellVertices,
    const AllCellFaces& cellFaces,
    const AllCellEdges& cellEdges,
    std::vector<spCell>& gCells
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
    const AllCellVertices& cellVertices,
    const AllCellFaces& cellFaces,
    const AllCellEdges& cellEdges,
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
        if (cell->convexes.size() > 0) {
            gCompounds.push_back(Compound{ cell->convexes });
        }
    }
}

////////////////////////////////////////////////////////////////////
//TEST for Customize Readobj function
////////////////////////////////////////////////////////////////////
void testObj(const std::string& filePath, 
    std::vector<spConvex>& compound) {
    compound = readOBJByComponents(filePath);
}

////////////////////////////////////////////////////////////////////
//TEST for Pipeline 
////////////////////////////////////////////////////////////////////
void testPipeline(const std::string& filePath,
    Pattern& pattern, 
    std::vector<Compound>& splitted) {
    auto convexes = readOBJByComponents(filePath);

    // Using weighted sum to approximate compound's CoM
    // Some of the original convex hulls might be overlapping so this is just an approximate.
    Eigen::Vector3d centroid = calculateCentroidCompound(convexes);
    Compound original{convexes, centroid};
    
    splitted = fracturePipeline(original, pattern);
}