if(TARGET igl::core)
    return()
endif()

include(FetchContent)
FetchContent_Declare(
    libigl
    GIT_REPOSITORY https://github.com/libigl/libigl.git
    GIT_TAG c7a84522c328a8e396205bbe92c40f727d709c8f # latest commit by 3/18/2024
)
FetchContent_MakeAvailable(libigl)
