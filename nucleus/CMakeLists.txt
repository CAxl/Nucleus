cmake_minimum_required(VERSION 3.21)

project(nucleus)

# Collect all cpp files under source/
file(GLOB_RECURSE SourceFiles "source/*.cpp")

# Build the executable
add_executable(${PROJECT_NAME} ${SourceFiles})

# Tell the compiler where to find headers
target_include_directories(${PROJECT_NAME}
    PUBLIC
        ${PROJECT_SOURCE_DIR}/include       # for <shells/...>, <observables/...>
        ${PROJECT_SOURCE_DIR}/source        # so "PhysConstants.h" still works
)

# Link against your math library
target_link_libraries(${PROJECT_NAME} PRIVATE math)
