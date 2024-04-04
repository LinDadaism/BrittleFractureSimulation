include(FetchContent)
FetchContent_Declare(
    geogram
    GIT_REPOSITORY https://github.com/BrunoLevy/geogram.git
    GIT_TAG v1.8.8
    OVERRIDE_FIND_PACKAGE
)
FetchContent_MakeAvailable(geogram)
