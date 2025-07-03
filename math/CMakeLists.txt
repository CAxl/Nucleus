cmake_minimum_required(VERSION 3.21)

project(math)
file(GLOB_RECURSE SourceFiles "source/*.cpp")

add_library("${PROJECT_NAME}" "${SourceFiles}")



set(EigenLib "${PROJECT_SOURCE_DIR}/libs/")

target_include_directories("${PROJECT_NAME}" PRIVATE "${PROJECT_SOURCE_DIR}/source", 
PUBLIC "${PROJECT_SOURCE_DIR}/include" 
"${EigenLib}"

)

