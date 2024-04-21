#pragma once
#include <fstream>
#include <sstream>
#include <string>
#include <stdexcept>
#include "clipper.h"

// Hard-coded 10 Voronoi cells as fracture pattern for Maya Plugin
const static AllCellVertices gCellVertices10   // Vertices of each Voronoi cell
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
const static AllCellFaces gCellFaces10           // Faces of each Voronoi cell represented by vertex indices
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
const static AllCellEdges gCellEdges10              // Edges of each Voronoi cell represented by vertex indices
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

std::vector<spConvex> readOBJByComponents(const std::string& filePath);
