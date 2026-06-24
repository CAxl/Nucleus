cmake_minimum_required(VERSION 3.21)

project(math)

file(GLOB_RECURSE SourceFiles "source/*.cpp")
add_library(${PROJECT_NAME} ${SourceFiles})

target_include_directories(${PROJECT_NAME}
    PUBLIC
        ${PROJECT_SOURCE_DIR}/include
        ${PROJECT_SOURCE_DIR}/libs
        ${PROJECT_SOURCE_DIR}/libs/Spectra/include
    PRIVATE
        ${PROJECT_SOURCE_DIR}/source
)
