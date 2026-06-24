cmake_minimum_required(VERSION 3.21)

project(nucleus)

file(GLOB_RECURSE SourceFiles "source/*.cpp")
add_executable(${PROJECT_NAME} ${SourceFiles})

target_link_libraries(${PROJECT_NAME} PRIVATE math)

target_include_directories(${PROJECT_NAME}
    PUBLIC
        ${PROJECT_SOURCE_DIR}/include
)
