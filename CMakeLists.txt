cmake_minimum_required(VERSION 3.16)
project(fracture)

list(PREPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

# ============================================================================== 
# Libigl
# ==============================================================================
include(libigl)

# Enable the target igl::glfw
igl_include(glfw)
# Other modules you could enable
#igl_include(embree)
igl_include(imgui)
igl_include(opengl)
#igl_include(stb)
#igl_include(predicates)
#igl_include(xml)
igl_include(copyleft cgal)
#igl_include(copyleft comiso)
#igl_include(copyleft core)
#igl_include(copyleft cork)
igl_include(copyleft tetgen)
#igl_include(restricted matlab)
#igl_include(restricted mosek)
#igl_include(restricted triangle)

# Add your project files
file(GLOB SRC_FILES src/main.cpp
		    src/vec.cpp src/vec.h
	 	    src/clipper.cpp src/clipper.h)
add_executable(${PROJECT_NAME} ${SRC_FILES})

# Link igl (and the glfw module) to your project
target_link_libraries(${PROJECT_NAME} PUBLIC 
  igl::glfw
  ## Other modules you could link to
  # igl::embree
  igl::imgui
  igl::opengl
  # igl::stb
  # igl::predicates
  # igl::xml
  igl_copyleft::cgal
  # igl_copyleft::comiso
  # igl_copyleft::core
  # igl_copyleft::cork
  igl_copyleft::tetgen
  # igl_restricted::matlab
  # igl_restricted::mosek
  # igl_restricted::triangle
)

# ==============================================================================
# Voro++
# ==============================================================================
# List of headers locations
include_directories("voro/include")

# List of cpu sources
file(GLOB_RECURSE host_sources voro/src/*.cc)

add_library(voro++ STATIC ${host_sources})

target_link_libraries(${PROJECT_NAME} PRIVATE voro++)


# ==============================================================================
# GEOGRAM
# ==============================================================================
# set(CMAKE_TOOLCHAIN_FILE "C:/dev/vcpkg/scripts/buildsystems/vcpkg.cmake"
#  CACHE STRING "Vcpkg toolchain file")
# include(C:/dev/vcpkg/scripts/buildsystems/vcpkg.cmake)
# find_package(Geogram REQUIRED)
# message(STATUS "Geogram: ${GEOGRAM_LIBRARY}, ${GEOGRAM_INCLUDE_DIR}")
# target_link_libraries(${PROJECT_NAME} PRIVATE Geogram::geogram)

# no need to install manually
# include(geogram)
# target_link_libraries(${PROJECT_NAME} PRIVATE geogram)