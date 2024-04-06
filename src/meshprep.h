#ifndef MESH_PREP
#define MESH_PREP
#include <fstream>
#include <sstream>
#include <string>
#include <stdexcept>
#include "clipper.h"

Compound readOBJByComponents(const std::string& filePath);


#endif