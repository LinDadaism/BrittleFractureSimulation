
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
* [Voro++](https://math.lbl.gov/voro++/), a C++ library for carrying out three-dimensional computations of the Voronoi tessellation.
* [CGAL](https://www.cgal.org/), a C++ library of computational geometry algorithms mainly for mesh clipping.
* [CoACD](https://github.com/SarahWeiii/CoACD), an approximate convex decomposition algorithm that we use for mesh preparation before fracture.

**We configure the project using cmake. Installation of `libigl` is handled automatically. We include source code of `voro++` in the project, thus building `voro++` also automatically handled by cmake.** 

Now let's build the Visual Studio project:

1. From the `Brittle-Fracture` directory, run:

    * Windows: `cmake -B build -S . -G "Visual Studio 17 2022" -A x64` to generate the Visual Studio project. You may choose a different Visual Studio version.
2. Open the generated `BFX.sln` project inside `/build` directory.
3. Build and Run. Make sure you run the `BFXViewer` target by right-clicking it and selecting "Set as StartUp Project".

A GUI window should launch displaying a 3D cube along with some debugging visuals of the vornoi decomposition.


### Troubleshooting
When you build the project for the **first** time (or re-build), there's an error:
```
...\BFX\voro\src\v_base_wl.cc(17,65): error C2888: 'const unsigned int voro::voro_base::wl[4096]': symbol cannot be defined within namespace 'voro'
```
This is a weird issue in the `voro++` library, and the trick we found to bypass the error is also very weird. Do the following (order matters):

1. Comment out the 3 lines below in the file and build the project (don't rebuild).
```
#include "v_base.hh"
namespace voro {
    ...
}
```
2. Uncomment the header line `#include "v_base.hh"`, and build the project (don't rebuild).
3. Uncomment the namespace parentheses `namespace voro{ ... }`, and build the project (don't rebuild).

There should be no more errors. `voro++.lib` and `voro++.vcxproj` should be successfully built.

