/**
 * This file is more of math utility functions.
**/
#pragma once

#include <Windows.h>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include "vec.h"

/////////////////////////////////////////////////////////////
///             Matrix - Array Conversions                 //
/////////////////////////////////////////////////////////////

Eigen::MatrixXd convertToMatrixXd2D(const std::vector<std::vector<Eigen::Vector3d>>& vectorOfVectors) {
    // First, calculate the total number of rows required in the matrix
    int totalRows = 0;
    for (const auto& vec : vectorOfVectors) {
        totalRows += vec.size();
    }

    // Create an Eigen::MatrixXd with the correct size
    Eigen::MatrixXd result(totalRows, 3); // 3 columns for x, y, z coordinates

    // Fill the matrix
    int currentRow = 0;
    for (const auto& vec : vectorOfVectors) {
        for (const auto& point : vec) {
            result(currentRow, 0) = point.x();
            result(currentRow, 1) = point.y();
            result(currentRow, 2) = point.z();
            currentRow++;
        }
    }

    return result;
}

Eigen::MatrixXd convertToMatrixXd(const std::vector<Eigen::Vector3d>& vec) {
    Eigen::MatrixXd mat(vec.size(), 3); // 3 columns for x, y, z
    for (int i = 0; i < vec.size(); ++i) {
        mat.row(i) = vec[i];
    }

    return mat;
}

Eigen::MatrixXi convertToMatrixXi(const std::vector<std::vector<int>>& vec) {
    Eigen::MatrixXi mat(vec.size(), 3); // 3 columns for x, y, z
    for (int i = 0; i < vec.size(); ++i) {
        mat(i, 0) = vec[i][0];
        mat(i, 1) = vec[i][1];
        mat(i, 2) = vec[i][2];
    }

    return mat;
}

void convertToFaceVertArray(const Eigen::MatrixXi& faceVertsMat, std::vector<std::vector<int>>& faceVerts)
{
    for (int i = 0; i < faceVertsMat.rows(); ++i) {
        // Temporary vector to hold the vertex indices for the current face
        std::vector<int> indices;

        for (int j = 0; j < faceVertsMat.cols(); ++j) {
            indices.push_back(faceVertsMat(i, j));
        }
        faceVerts.push_back(indices);
    }
}

void convertToVertArray(const Eigen::MatrixXd& vertMatrix, std::vector<Eigen::Vector3d>& vertVector)
{
    for (int i = 0; i < vertMatrix.rows(); ++i) {
        Eigen::Vector3d vertex = vertMatrix.row(i);
        vertVector.push_back(vertex);
    }
}

void convertEigenToVec(const std::vector<Eigen::Vector3d>& eigenArr, std::vector<vec3>& vecArr) {
    vecArr.clear();
    for (const auto& eigenVec : eigenArr) {
        vecArr.push_back(vec3(eigenVec.x(), eigenVec.y(), eigenVec.z()));
    }
}

// decompose an AABB into 8 smaller, equally sized cells
void decomposeAABB(const Eigen::Vector3d& minCorner, 
    const Eigen::Vector3d& maxCorner,
    std::vector<std::vector<Eigen::Vector3d>>& cellVertices,
    std::vector<std::vector<std::vector<int>>>& cellFaces
) {
    // Compute cell dimensions
    Eigen::Vector3d cellSize = (maxCorner - minCorner) / 2.0;

    // Define faces for each cell using vertex indices
    std::vector<std::vector<int>> faces = {
        // Bottom face
        {0, 1, 2}, {0, 2, 3},
        // Top face
        {4, 5, 6}, {4, 6, 7},
        // Front face
        {0, 1, 5}, {0, 5, 4},
        // Back face
        {2, 3, 7}, {2, 7, 6},
        // Left face
        {0, 3, 7}, {0, 7, 4},
        // Right face
        {1, 2, 6}, {1, 6, 5}
    };

    // Compute vertices for each of the 8 cells
    for (int i = 0; i < 2; ++i) {
        for (int j = 0; j < 2; ++j) {
            for (int k = 0; k < 2; ++k) {
                Eigen::Vector3d baseCorner = minCorner + Eigen::Vector3d(i * cellSize.x(), j * cellSize.y(), k * cellSize.z());
                std::vector<Eigen::Vector3d> vertices = {
                    baseCorner,
                    baseCorner + Eigen::Vector3d(cellSize.x(), 0, 0),
                    baseCorner + Eigen::Vector3d(cellSize.x(), cellSize.y(), 0),
                    baseCorner + Eigen::Vector3d(0, cellSize.y(), 0),
                    baseCorner + Eigen::Vector3d(0, 0, cellSize.z()),
                    baseCorner + Eigen::Vector3d(cellSize.x(), 0, cellSize.z()),
                    baseCorner + Eigen::Vector3d(cellSize.x(), cellSize.y(), cellSize.z()),
                    baseCorner + Eigen::Vector3d(0, cellSize.y(), cellSize.z())
                };
                cellVertices.push_back(vertices);
                
                // All cells have the same faces, so we assign the same face indices to each cell
                cellFaces.push_back(faces);
            }
        }
    }
}

// We use CreateProcess() to run an external executable 
// https://learn.microsoft.com/en-us/windows/win32/procthread/creating-processes
bool executeCommand(LPCTSTR lpApplicationName, LPSTR lpCommand)
{
    // additional information
    STARTUPINFO si;
    PROCESS_INFORMATION pi;

    // set the size of the structures
    ZeroMemory(&si, sizeof(si));
    si.cb = sizeof(si);
    ZeroMemory(&pi, sizeof(pi));

    // start the program up
    if (!CreateProcess(lpApplicationName,   // the module to be executed, can be NULL. In that case, the module name must be the first white space–delimited token in the lpCommand string. 
            lpCommand,      // Command line
            NULL,           // Process handle not inheritable
            NULL,           // Thread handle not inheritable
            FALSE,          // Set handle inheritance to FALSE
            0,              // No creation flags
            NULL,           // Use parent's environment block
            NULL,           // Use parent's starting directory 
            &si,            // Pointer to STARTUPINFO structure
            &pi)            // Pointer to PROCESS_INFORMATION structure (removed extra parentheses)
        ) 
    {
        printf("CreateProcess failed (%d).\n", GetLastError());
        return false;
    }

    // Wait until child process exits.
    WaitForSingleObject(pi.hProcess, INFINITE);

    // Close process and thread handles. 
    CloseHandle(pi.hProcess);
    CloseHandle(pi.hThread);

    return true;
}

void printPattern(vector<Eigen::Vector3d>& gPoints,
    std::vector<std::vector<Eigen::Vector3d>>& gCellVertices,
    std::vector<std::vector<std::vector<int>>>& gCellFaces,
    std::vector<Eigen::MatrixXi>& gCellEdges)
{
    // gPoints
    std::cout << "cellNodes:\n{" << std::endl;
    for (const auto& p : gPoints) {
        std::cout << "  Eigen::Vector3d(" << p.x() << ", " << p.y() << ", " << p.z() << ")," << std::endl;
    }
    std::cout << "}\n" << std::endl;

    // gCellVertices
    std::cout << "cellVerts:\n{" << std::endl;
    for (const auto& cell : gCellVertices) {
        std::cout << "  {" << std::endl;
        for (const auto& vert : cell) {
            std::cout << "      Eigen::Vector3d(" << vert.x() << ", " << vert.y() << ", " << vert.z() << ")," << std::endl;
        }
        std::cout << "  }," << std::endl;
    }
    std::cout << "}\n" << std::endl;

    // gCellFaces
    std::cout << "cellFaces:\n{" << std::endl;
    for (const auto& face : gCellFaces) {
        std::cout << "  {" << std::endl;
        for (const auto& vertIds : face) {
            std::cout << "      {";
            for (int i = 0; i < vertIds.size(); i++) {
                if (i == vertIds.size() - 1) {
                    std::cout << vertIds[i];
                    break;
                }
                std::cout << vertIds[i] << ", ";
            }
            std::cout << "}," << std::endl;
        }
        std::cout << "  }," << std::endl;
    }
    std::cout << "}\n" << std::endl;

    // gCellEdges
    std::cout << "cellEdges:\n{" << std::endl;
    for (const auto& cell : gCellEdges) {
        int numRows = static_cast<int>(cell.rows());
        std::cout << "  (Eigen::MatrixXi(" << numRows << ", 2) << ";
        for (int i = 0; i < cell.size(); i++) {
            if (i == cell.size() - 1) {
                std::cout << cell(i) << ").finished()," << std::endl;
                break;
            }
            std::cout << cell(i) << ", ";
        }
    }
    std::cout << "}\n" << std::endl;
}