
# BFX for VFX: Brittle Fracture Simulation

Authors: Yilin Li and Linda Zhu

## Installation

### Dependencies
BFX depends on external libraries:

* [libigl](https://libigl.github.io/), a simple geometry processing library used for mesh prep and visualization. We also use the following modules linked with `libigl`: 
    - imgui
    - opengl
    - glfw
    - tetgen
* [Geogram](https://github.com/BrunoLevy/geogram?tab=readme-ov-file), a geometry processing library mainly used for Voronoi decomposition and mesh intersection.

**We configure the project using cmake. Installation of `libigl` is handled automatically. However, we do need to manually install `geogram` to the local machine.** The following instructions walk through the recommended steps for installing `geogram`on Windows.

1. Install `vcpkg`

[Vcpkg](https://github.com/alicevision/vcpkg) is a package manager that helps in acquiring, building, and managing C/C++ libraries. `geogram` can be built with it. See the reference [installation guide](https://github.com/alicevision/vcpkg/blob/alicevision_master/README.md#quick-start-windows) to setup vcpkg.

```bash
cd C:\dev\ # (recommended but optional, reason being that 1. global install - somewhere like C:\src\vcpkg or C:\dev\vcpkg, avoids path issues for some port build systems; 2. this is the default path CMakeLists searches for the package)
git clone https://github.com/microsoft/vcpkg
.\vcpkg\bootstrap-vcpkg.bat
```
**Note**: If you downloaded `vcpkg` somewhere else, you have to change the `CMAKE_TOOLCHAIN_FILE` in the CMakeLists.txt to reflect that change!

2. Build the required dependencies
```bash
cd <VCPKG_INSTALL_DIR>
set VCPKG_ROOT=%cd%

vcpkg install geogram # This might take a while.

# geogram should show up if you call `vcpkg list`
```

Now we can build the Visual Studio project:

3. From the `Brittle-Fracture` directory, run:

    * Windows: `cmake -B build -S . -G "Visual Studio 17 2022" -A x64` to generate the Visual Studio project. You may choose a different Visual Studio version.
4. Open the generated `fracture.sln` project inside `/build` directory.
5. Build and Run. Make sure you run the `transpose` target (not `ALL_BUILD`) by right-clicking it and selecting "Set as StartUp Project".
6. If there's a compile error in `tetgenio_to_tetmesh.cpp`,
```
'numberofpointmarkers' is not a member of 'tetgenio'
```
**delete line 86 of this file (it's a bug)!**
```
assert(out.numberofpoints == out.numberofpointmarkers);
```

A GUI window should launch displaying a 3D cube.


# Misc below:
## Using other modules of libigl

This example project uses the `igl::opengl::glfw::Viewer`, therefore it requires
the glfw module of libigl. This shows up in the CMakeLists.txt 

```cmake
igl_include(glfw)
…
target_link_libraries(${PROJECT_NAME} PUBLIC igl::glfw)
```

Suppose you also wanted to use the triangle module in libigl. Then you would
change these to

```cmake
igl_include(glfw)
igl_include(restricted triangle)
…
target_link_libraries(${PROJECT_NAME} PUBLIC igl::glfw igl_restricted::triangle)
```

The "restricted" appears in this case because the triangle library has a more
restricted license than libigl. See other examples commented out in
CMakeLists.txt.


## Dependencies

The only dependencies are STL, Eigen, [libigl](http://libigl.github.io/libigl/) and the dependencies
of the `igl::opengl::glfw::Viewer` (OpenGL, glad and GLFW).

The CMake build system will automatically download libigl and its dependencies using
[CMake FetchContent](https://cmake.org/cmake/help/latest/module/FetchContent.html),
thus requiring no setup on your part.

### Use a local copy of libigl
You can use the CMake cache variable `FETCHCONTENT_SOURCE_DIR_LIBIGL` when configuring your CMake project for
the first time to aim it at a local copy of libigl instead.
```
cmake -DFETCHCONTENT_SOURCE_DIR_LIBIGL=<path-to-libigl> ..
```
When changing this value, do not forget to clear your `CMakeCache.txt`, or to update the cache variable
via `cmake-gui` or `ccmake`.
