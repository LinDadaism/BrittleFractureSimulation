#pragma once

#include "meshprep.h"

class Simulation
{
public:
    Simulation();
    Simulation(std::string filepath);
    Simulation(vec3 pos, float r, float amt, std::string filepath);
    ~Simulation() {};

    std::string getInputFilepath();
    std::vector<Compound> getFractureShards();
    int numPatternCells();

    void setInputFilepath(std::string filepath);

    // generate fractures based on hard-coded Voronoi patterns
    void genFractureUniformStatic();
    // generate fractures based on dynamically placed Voronoi cells
    void genFractureUniformDynamic(const std::vector<vec3>& nodes, vec3 minCorner, vec3 maxCorner);

private:
    // Fracture simulation config
    vec3 impactPos;
    float impactRadius = -1.f;      // default -1 means full-body fracture
    float explodeAmt = 1.f;         // distance between shards

    // Mesh info
    //vec3 meshPos;
    std::string coacdMeshFilepath;
    std::vector<Compound> fractureShards;

    // Pattern info
    Pattern pattern;
};