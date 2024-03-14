cmake_minimum_required(VERSION 3.10)

# set the project name
project(BFX)

# specify the C++ standard
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED True)

set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR})
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR})
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR})

# Set up libraries

add_library(particle STATIC
    ./src/animation/aParticle.h
    ./src/animation/aParticle.cpp
    ./src/animation/aParticleSystem.h
    ./src/animation/aParticleSystem.cpp
)

target_include_directories(particle PUBLIC
    ./3rdparty/eigen-3.2.8
    ./3rdparty/glm
    ./src/animation
)

# Set up executables/viewers
# Find OpenGL
find_package(OpenGL REQUIRED)

# Add GLFW
set(GLFW_BUILD_DOCS OFF CACHE BOOL "" FORCE)
set(GLFW_BUILD_TESTS OFF CACHE BOOL "" FORCE)
set(GLFW_BUILD_EXAMPLES OFF CACHE BOOL "" FORCE)
add_subdirectory(./3rdparty/glfw-3.3.2)

# Add GLAD
add_subdirectory(./3rdparty/glad)

# Add ImGUI
add_subdirectory(./3rdparty/imgui-1.77)

# Add tinyobj
add_library(tinyobj STATIC
    ./3rdparty/tinyobj/tiny_obj_loader.h
    ./3rdparty/tinyobj/tiny_obj_loader.cc
)

add_executable(ParticleViewer
    ./src/viewer/particleMain.cpp
    ./src/viewer/particleViewer.h
    ./src/viewer/particleViewer.cpp
    ./src/viewer/camera.cpp
    ./src/viewer/camera.h
    ./src/viewer/objmodel.h
    ./src/viewer/objmodel.cpp
    ./src/viewer/shader.h
    ./src/viewer/viewer.cpp
    ./src/viewer/viewer.h
    ./src/viewer/drawable.h
    ./src/viewer/drawable.cpp
)

target_include_directories(ParticleViewer PUBLIC 
    ./src/viewer
    ./src/animation
    ./3rdparty/imgui-1.77
    ./3rdparty/glm
    ./3rdparty/tinyobj
)

target_link_libraries(ParticleViewer PUBLIC particle glad glfw imgui tinyobj)