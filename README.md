
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
    - cgal
* [Voro++](https://math.lbl.gov/voro++/), a C++ library for carrying out three-dimensional computations of the Voronoi tessellation.

**We configure the project using cmake. Installation of `libigl` is handled automatically. We include source code of `voro++` in the project, thus building `voro++` also automatically handled by cmake.** 

Now let's build the Visual Studio project:

3. From the `Brittle-Fracture` directory, run:

    * Windows: `cmake -B build -S . -G "Visual Studio 17 2022" -A x64` to generate the Visual Studio project. You may choose a different Visual Studio version.
4. Open the generated `fracture.sln` project inside `/build` directory.
5. Build and Run. Make sure you run the `transpose` target (not `ALL_BUILD`) by right-clicking it and selecting "Set as StartUp Project".
6. **ATTENTION!** If you are using the `CGAL` module wrapper of `libigl`, you need to move the files in `/gmp` into the following build directories:

    - Cut/paste the `/gmp-src` folder into `/build/_deps`
    - Cut/paste the `/mpfr-src` folder into `/build/_deps`

A GUI window should launch displaying a 3D cube along with some debugging visuals of the vornoi decomposition.
