cmake_minimum_required(VERSION 3.16)
project(OrbitalVisualizer)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# ── Find system libraries (installed via apt) ────────────────────────────────
find_package(OpenGL REQUIRED)
find_package(GLEW   REQUIRED)
find_package(glfw3  REQUIRED)

# GLM is header-only; find_package gives us the include path
find_package(glm    QUIET)
if(NOT glm_FOUND)
    # Fallback: glm headers are usually at /usr/include/glm on Ubuntu
    find_path(GLM_INCLUDE_DIR glm/glm.hpp REQUIRED)
    add_library(glm::glm INTERFACE IMPORTED)
    set_target_properties(glm::glm PROPERTIES
        INTERFACE_INCLUDE_DIRECTORIES "${GLM_INCLUDE_DIR}")
endif()

# ── Helper to register one visualizer executable ─────────────────────────────
function(add_visualizer target source)
    add_executable(${target} src/${source})
    target_link_libraries(${target} PRIVATE
        OpenGL::GL
        GLEW::GLEW
        glfw
        glm::glm
    )
    # Some GL headers need this on Linux
    target_compile_definitions(${target} PRIVATE GL_SILENCE_DEPRECATION)
endfunction()

# ── Targets ──────────────────────────────────────────────────────────────────
add_visualizer(atom            atom.cpp)
add_visualizer(atom_realtime   atom_realtime.cpp)
add_visualizer(atom_raytracer  atom_raytracer.cpp)
add_visualizer(wave_atom_2d    wave_atom_2d.cpp)
