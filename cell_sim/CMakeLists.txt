cmake_minimum_required(VERSION 3.10)  # CMake version check
project(capmd)                        # Create project "capmd"
set(CMAKE_CXX_STANDARD 14)            # Enable c++14 standard
set(CMAKE_CXX_STANDARD_REQUIRED ON)   # Require c++14 standard
#[[
if(APPLE AND EXISTS /usr/local/opt/qt5)
    # Homebrew installs Qt5 (up to at least 5.9.1) in /usr/local/qt5.
    # This check makes sure that it's located rather than the default /usr/local prefix.
    list(APPEND CMAKE_PREFIX_PATH "/usr/local/opt/qt5")
endif()
if(APPLE AND EXISTS /usr/local/Cellar/vtk@8.2/8.2.0_6)
    # Even worse, homebrew installs new versions of vtk which are not longer supported
    # Instead, use an old version and make explicit here
    set(VTK_DIR "/usr/local/Cellar/vtk@8.2/8.2.0_6/lib/cmake/vtk-8.2" CACHE PATH "VTK directory override" FORCE)
    message(STATUS "include dir: ${VTK_DIR}")
endif()
if (WIN32)
    # Same problem of versions for windows builds
    set(VTK_DIR "C:/MyProjects/VTK/lib/cmake/vtk-8.2" CACHE PATH "VTK directory override" FORCE)
endif (WIN32)

find_package(VTK REQUIRED NO_MODULE)
include(${VTK_USE_FILE})
]]#

find_package(PythonLibs REQUIRED)  # Require Python librarys
include_directories(${PYTHON_INCLUDE_DIRS})  # include python directories
add_library(
    capmd
    ${CMAKE_SOURCE_DIR}/src/system/Simulation.cpp
    ${CMAKE_SOURCE_DIR}/src/system/Particle.cpp
    ${CMAKE_SOURCE_DIR}/src/system/Population.cpp
    ${CMAKE_SOURCE_DIR}/src/system/Domain.cpp
    ${CMAKE_SOURCE_DIR}/src/system/Dynamics.cpp
    ${CMAKE_SOURCE_DIR}/src/system/potentials/Interaction.cpp
    ${CMAKE_SOURCE_DIR}/src/system/potentials/Potential.cpp
    ${CMAKE_SOURCE_DIR}/src/system/potentials/newPotential.cpp
    ${CMAKE_SOURCE_DIR}/src/system/potentials/adhPotential.cpp
    ${CMAKE_SOURCE_DIR}/src/bindings/Interface.cpp
    ${CMAKE_SOURCE_DIR}/src/io/Output.cpp)

# Add additonal directories.
include_directories(include/config)
include_directories(src/system)
include_directories(src/system/potentials)
include_directories(src/io)
include_directories(src/bindings)

# Add the testing subdir.
add_subdirectory(tests)

target_compile_options(capmd PRIVATE -fpie)
#target_compile_options(capmd PRIVATE -fPIC)

# Add executable target for library demo.
add_executable(
    capmd_demo
    ${CMAKE_SOURCE_DIR}/src/main.cpp
)

# Link the CAPMD demo against libcapmd, VTK libraries
target_link_libraries(
        capmd_demo PRIVATE
        capmd
#        ${VTK_LIBRARIES}
)

# Add pybind11 folder and make python module (called pycapmd)
add_subdirectory(libs/pybind11 ${CMAKE_SOURCE_DIR}/libs/pybind11)
pybind11_add_module(
    pycapmd
    ${CMAKE_SOURCE_DIR}/src/bindings/bindings.cpp
)
set_target_properties(pycapmd PROPERTIES
        LIBRARY_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/../simulator/)

# Link vtk libraries to python module
target_link_libraries(
        pycapmd PUBLIC
        capmd
        #${VTK_LIBRARIES}
)
