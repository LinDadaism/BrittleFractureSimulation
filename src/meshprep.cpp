#include "meshprep.h"

Compound readOBJByComponents(const std::string& filePath) {
    std::ifstream file(filePath);
    std::vector<spConvex> meshes;
    std::string line;

    if (!file.is_open()) {
        throw std::runtime_error("Could not open file: " + filePath);
    }
    
    int currentVertices = 0; 
    std::vector<Eigen::Vector3d> vertices; 
    std::vector<std::vector<int>> faces;
    while (std::getline(file, line)) {
        std::istringstream iss(line);
        std::string type;
        iss >> type;
        if (type == "o" || type == "g") {
            // Start a new mesh on encountering a new object/group name
            if (vertices.size() > 0) {
                Surface_mesh sm;
                buildSMfromVF(vertices, faces, sm);
                spConvex cur(new MeshConvex{ vertices, faces, sm});
                
                // Update centroid for the current mesh
                calculateCentroid(*cur, Eigen::Vector3d(0, 0, 0));

                meshes.push_back(cur);
                currentVertices += vertices.size();
                
                vertices.clear();
                faces.clear();
            }
        }
        else if (type == "v") {
            double x, y, z;
            iss >> x >> y >> z;
            vertices.push_back(Eigen::Vector3d(x, y, z));
        }
        else if (type == "f") {
            int a, b, c;
            iss >> a >> b >> c; 
            std::vector<int> face{ a - 1 - currentVertices, b - 1 - currentVertices, c - 1 - currentVertices };
            faces.push_back(face);
        }
    }
    if (vertices.size() > 0) {
        Surface_mesh sm;
        buildSMfromVF(vertices, faces, sm);
        spConvex cur(new MeshConvex{ vertices, faces, sm });
        // Update centroid for the current mesh
        calculateCentroid(*cur, Eigen::Vector3d(0, 0, 0));
        meshes.push_back(cur);
    }
    return Compound{ meshes };
}
