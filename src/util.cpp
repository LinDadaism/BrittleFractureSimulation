#include "util.h"

void meshClipping(Plane plane, const Eigen::MatrixXd meshV, const Eigen::MatrixXi meshF) {
//    std::vector<Eigen::Vector3d> newVertices;
//    std::vector<std::vector<int>> newFaces;
//
//    for (unsigned int i = 0; i < meshF.rows(); ++i) {
//        std::vector<int> newFace;
//        for (unsigned int j = 0; j < meshF.cols(); ++j) {
//            int current = meshF[i, j];
//            int next = meshF[i, (j+1) % meshF.cols()];
//
//            vec3 curVertice = vec3(meshV[current, 0], meshV[current, 1], meshV[current, 2]);
//            vec3 nextVertice = vec3(meshV[next, 0], meshV[next, 1], meshV[next, 2]);
//            float currentDistance = plane.distanceToPoint(curVertice);
//            float nextDistance = plane.distanceToPoint(nextVertice);
//
//            if (currentDistance >= 0) {
//                newFaces.push_back(mesh.addVertex(currentVertex)); // Keep the original vertex if it's on the visible side of the plane
//            }
//
//            if (currentDistance * nextDistance < 0) { // Edge crosses the plane
//                // Compute intersection point
//                float t = currentDistance / (currentDistance - nextDistance);
//                Vector3 intersection = currentVertex + (nextVertex - currentVertex) * t;
//
//                // Add intersection point as a new vertex
//                newFace.push_back(mesh.addVertex(intersection));
//            }
//        }
//
//        if (!newFace.empty()) {
//            newFaces.push_back(newFace);
//        }
//    }
//
//    // Update the mesh with the new vertices and faces
//    mesh.faces = newFaces;
}
