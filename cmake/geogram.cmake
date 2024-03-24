include(FetchContent)
FetchContent_Declare(
    geogram
    GIT_REPOSITORY https://github.com/BrunoLevy/geogram.git
    GIT_TAG v1.8.8
)
FetchContent_MakeAvailable(geogram)
