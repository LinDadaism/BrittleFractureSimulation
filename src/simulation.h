#pragma once

#include "meshprep.h"

class Simulation
{
public:
    Simulation();
    Simulation(std::string filepath);
    Simulation(Eigen::Vector3d pos, float r, float amt, std::string filepath);
    ~Simulation() {};

    std::string getInputFilepath();
    std::vector<Compound> getFractureShards();
    int numPatternCells();

    void setInputFilepath(std::string filepath);

    // generate fractures based on hard-coded Voronoi patterns
    void genFractureUniformStatic();
    // generate fractures based on dynamically placed Voronoi cells
    void genFractureUniformDynamic(const std::vector<Eigen::Vector3d>& nodes, Eigen::Vector3d minCorner, Eigen::Vector3d maxCorner);

private:
    // Fracture simulation config
    Eigen::Vector3d impactPos;
    float impactRadius = -1.f;      // default -1 means full-body fracture
    float explodeAmt = 1.f;         // distance between shards

    // Mesh info
    //vec3 meshPos;
    std::string coacdMeshFilepath;
    std::vector<Compound> fractureShards;

    // Pattern info
    Pattern pattern;
};