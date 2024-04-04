/**
 * This file is more of math utility functions.
**/
#pragma once

#include <Eigen/Dense>
#include <Eigen/Sparse>

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