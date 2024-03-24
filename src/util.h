#ifndef UTIL_FILE
#define UTIL_FILE
#include "vec.h"
#include <Eigen/Dense>
#include <Eigen/Sparse>

struct Plane {
    vec3 normal;
    float distance;

    Plane(const vec3& normal, float distance) : normal(normal), distance(distance) {}

    // Function to compute distance from the plane to a point
    float distanceToPoint(const vec3& point) const {
        return normal[0] * point[0] + normal[1] * point[1] + normal[2] * point[2] - distance;
    }
};

void meshClipping(Plane, const Eigen::MatrixXd, const Eigen::MatrixXi, std::vector<Eigen::Vector3d>, vector<vector<int>>);

#endif // UTIL_FILE
