set(CMAKE_BUILD_TYPE "Release" CACHE INTERNAL "Build a release version of CoACD" FORCE)
set(CMAKE_MSVC_RUNTIME_LIBRARY "MultiThreaded" CACHE INTERNAL "" FORCE)
set(OPENVDB_CORE_SHARED OFF CACHE BOOL "" FORCE)
set(TBB_TEST OFF CACHE BOOL "" FORCE)
set(CMAKE_CXX_FLAGS "/MT /EHsc" CACHE STRING "" FORCE)
set(WITH_3RD_PARTY_LIBS OFF CACHE BOOL "Disable 3rd party libs in CoACD" FORCE)

include(FetchContent)
FetchContent_Declare(
    coacd
    GIT_REPOSITORY https://github.com/SarahWeiii/CoACD.git
    GIT_TAG 1.0.1
    GIT_PROGRESS TRUE
)

FetchContent_MakeAvailable(coacd)
