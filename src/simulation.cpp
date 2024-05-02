#include "simulation.h"

Simulation::Simulation()
{
}

Simulation::Simulation(std::string filepath) :
    coacdMeshFilepath(filepath)
{
}

Simulation::Simulation(Eigen::Vector3d pos, float r, float amt, std::string filepath) :
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

int Simulation::numPatternCells()
{
    return pattern.numCells();
}

void Simulation::setInputFilepath(std::string filepath)
{
	coacdMeshFilepath = filepath;
}

void Simulation::genFractureUniformStatic()
{
    pattern.setVertices(gCellVertices10);
    pattern.setFaces(gCellFaces10);
    pattern.setMin(Eigen::Vector3d(-0.5, -0.5, -0.5)); // TODO: update if pre-defined pattern changes
    pattern.setMax(Eigen::Vector3d(0.5, 0.5, 0.5));
    pattern.createCellsfromVoro();
    
    auto convexes = readOBJByComponents(coacdMeshFilepath);
   
    Eigen::Vector3d centroid = calculateCentroidCompound(convexes);
    Compound original{ convexes, centroid };

    fractureShards = fracturePipeline(original, pattern);
}

void Simulation::genFractureUniformDynamic(const std::vector<Eigen::Vector3d>& nodes, 
    Eigen::Vector3d minCorner, Eigen::Vector3d maxCorner)
{
    AllCellVertices cellVertices;   // Vertices of each Voronoi cell
    AllCellFaces cellFaces;         // Faces of each Voronoi cell represented by vertex indices
    AllCellEdges cellEdges;         // Edges of each Voronoi cell represented by vertex indices
    computeVoronoiCells(nodes, minCorner, maxCorner, cellVertices, cellFaces, cellEdges);

    pattern.setVertices(cellVertices);
    pattern.setFaces(cellFaces);
    pattern.setMin(minCorner);
    pattern.setMax(maxCorner);

    pattern.createCellsfromVoro();

    auto convexes = readOBJByComponents(coacdMeshFilepath);

    // Using weighted sum to approximate compound's CoM
    // Some of the original convex hulls might be overlapping so this is just an approximate.
    Eigen::Vector3d centroid = calculateCentroidCompound(convexes);
    Compound original{ convexes, centroid };

    fractureShards = fracturePipeline(original, pattern);
}