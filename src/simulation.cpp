#include "simulation.h"

Simulation::Simulation()
{
}

Simulation::Simulation(std::string filepath) : coacdMeshFilepath(filepath)
{
}

Simulation::Simulation(vec3 pos, float r, float amt, std::string filepath) :
	impactPos(pos), impactRadius(r), explodeAmt(amt),
	coacdMeshFilepath(filepath)
{
}

std::string Simulation::getInputFilepath()
{
	return coacdMeshFilepath;
}

std::vector<Compound> Simulation::getFractureShards()
{
	return fractureShards;
}

void Simulation::setInputFilepath(std::string filepath)
{
	coacdMeshFilepath = filepath;
}

void Simulation::genFractureUniform(const std::vector<vec3>& nodes, vec3 minCorner, vec3 maxCorner)
{
    AllCellVertices cellVertices;   // Vertices of each Voronoi cell
    AllCellFaces cellFaces;         // Faces of each Voronoi cell represented by vertex indices
    AllCellEdges cellEdges;         // Edges of each Voronoi cell represented by vertex indices
    computeVoronoiCells(nodes, minCorner, maxCorner, cellVertices, cellFaces, cellEdges);

    pattern.setVertices(cellVertices);
    pattern.setFaces(cellFaces);
    pattern.createCellsfromVoro();

    auto convexes = readOBJByComponents(coacdMeshFilepath);

    // Using weighted sum to approximate compound's CoM
    // Some of the original convex hulls might be overlapping so this is just an approximate.
    Eigen::Vector3d centroid = calculateCentroidCompound(convexes);
    Compound original{ convexes, centroid };

    fractureShards = fracturePipeline(original, pattern);
}