cmake_minimum_required(VERSION 3.21)

project(math)
file(GLOB_RECURSE SourceFiles "source/*.cpp")

add_library(${PROJECT_NAME} ${SourceFiles})

# Public include dirs: accessible to anything that links against math
target_include_directories(${PROJECT_NAME}
    PUBLIC
        ${PROJECT_SOURCE_DIR}/include
        ${PROJECT_SOURCE_DIR}/libs
        ${PROJECT_SOURCE_DIR}/libs/Spectra/include
)