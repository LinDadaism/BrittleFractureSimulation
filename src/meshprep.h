#ifndef MESH_PREP
#define MESH_PREP
#include <fstream>
#include <sstream>
#include <string>
#include <stdexcept>
#include "clipper.h"

std::vector<spConvex> readOBJByComponents(const std::string& filePath);


#endif