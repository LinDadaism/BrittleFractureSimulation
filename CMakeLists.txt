cmake_minimum_required(VERSION 3.16)
project(BFX)

list(PREPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin")
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib")
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib")

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
#igl_include(copyleft cgal)
#igl_include(copyleft comiso)
#igl_include(copyleft core)
#igl_include(copyleft cork)
igl_include(copyleft tetgen)
#igl_include(restricted matlab)
#igl_include(restricted mosek)
#igl_include(restricted triangle)

# Add your project files
file(GLOB_RECURSE BFX_SOURCE_SRCS ${CMAKE_CURRENT_SOURCE_DIR}/src/*)
list(REMOVE_ITEM BFX_SOURCE_SRCS ${CMAKE_CURRENT_SOURCE_DIR}/src/main.cpp)

file(GLOB SRC_FILES ${BFX_SOURCE_SRCS})
add_library(${PROJECT_NAME} STATIC ${SRC_FILES})

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
  # igl_copyleft::cgal
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
# include_directories("voro/include")

# # List of cpu sources
# file(GLOB_RECURSE host_sources voro/src/*.cc)

# add_library(voro++ STATIC ${host_sources})

# target_link_libraries(${PROJECT_NAME} PUBLIC voro++)

# ==============================================================================
# CGAL
# ==============================================================================
# specify toolchain location for cross-compiling
set(CMAKE_TOOLCHAIN_FILE "C:/dev/vcpkg/scripts/buildsystems/vcpkg.cmake" # TODO: make sure this matches the vcpkg path on your local machine
 CACHE STRING "Vcpkg toolchain file")
include(C:/dev/vcpkg/scripts/buildsystems/vcpkg.cmake) # TODO: make sure this matches the vcpkg path on your local machine

# List of headers locations
find_package(CGAL REQUIRED)
target_link_libraries(${PROJECT_NAME} PUBLIC CGAL::CGAL)

# Set executable test app
add_executable(BFXViewer src/main.cpp)

# Set additional configurations for the executable
target_include_directories(BFXViewer PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/src>
  $<INSTALL_INTERFACE:src>  # for clients using the installed library
)
target_link_libraries(BFXViewer PUBLIC BFX)