#ifndef MESH_PREP
#define MESH_PREP
#include <fstream>
#include <sstream>
#include <string>
#include <stdexcept>
#include "clipper.h"

namespace mesh_prep {
	std::vector<spConvex> readOBJ(const std::string& filePath);
}

#endif