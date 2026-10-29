cmake_minimum_required(VERSION 3.12...3.31)
project(Mesh_smoothing_3_Benchmark)

find_package(CGAL REQUIRED)
find_package(Eigen3 QUIET)
include(CGAL_Eigen3_support)

find_path(GMSH_INCLUDE_DIR NAMES gmsh.h
    HINTS ENV GMSH_ROOT "$ENV{HOME}/.local"
    PATH_SUFFIXES include)
if(GMSH_INCLUDE_DIR)
    # The Python SDK ships a versioned library without an unversioned symlink.
    file(STRINGS "${GMSH_INCLUDE_DIR}/gmsh.h" GMSH_VERSION_LINE
        REGEX "^#define GMSH_API_VERSION \"")
    string(REGEX MATCH "[0-9]+\\.[0-9]+" GMSH_API_VERSION "${GMSH_VERSION_LINE}")
    get_filename_component(GMSH_SDK_PREFIX "${GMSH_INCLUDE_DIR}" DIRECTORY)
endif()
find_library(GMSH_LIBRARY NAMES gmsh "libgmsh.so.${GMSH_API_VERSION}"
    NAMES_PER_DIR
    HINTS ENV GMSH_ROOT "${GMSH_SDK_PREFIX}" "$ENV{HOME}/.local"
    PATH_SUFFIXES lib lib64)
if(NOT GMSH_INCLUDE_DIR OR NOT GMSH_LIBRARY)
    message(FATAL_ERROR
        "Gmsh SDK not found. Set GMSH_INCLUDE_DIR and GMSH_LIBRARY, or add its prefix to CMAKE_PREFIX_PATH.")
endif()
add_library(Gmsh::Gmsh UNKNOWN IMPORTED)
set_target_properties(Gmsh::Gmsh PROPERTIES
    IMPORTED_LOCATION "${GMSH_LIBRARY}"
    INTERFACE_INCLUDE_DIRECTORIES "${GMSH_INCLUDE_DIR}")

add_library(benchmark_common
    common/cgal_common.cpp
    common/cgal_common.h
    common/distance_computation.cpp
    common/distance_computation.h
    common/gmsh_wrapper.cpp
    common/gmsh_wrapper.h
    common/statistics.h
    common/parameters.h
)
target_include_directories(benchmark_common PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/common")
target_link_libraries(benchmark_common PUBLIC CGAL::CGAL CGAL::Eigen3_support)
target_link_libraries(benchmark_common PRIVATE Gmsh::Gmsh)

find_package(TBB QUIET)
include (CGAL_TBB_support)
if(TARGET CGAL::TBB_support)
    target_link_libraries(benchmark_common PUBLIC CGAL::TBB_support)
endif()

create_single_source_cgal_program("common/result_header.cpp")
target_link_libraries(result_header PRIVATE benchmark_common)

create_single_source_cgal_program("SDF_meshing/sdf_meshing.cpp")
target_link_libraries(sdf_meshing PRIVATE benchmark_common)

create_single_source_cgal_program("octree_hex_meshing/octree_hex_mesh.cpp")
target_link_libraries(octree_hex_mesh PRIVATE benchmark_common)

create_single_source_cgal_program("offset_AW/offset_meshing.cpp")
target_link_libraries(offset_meshing PRIVATE benchmark_common)

create_single_source_cgal_program("poly_meshing/poly_remesh.cpp")
target_link_libraries(poly_remesh PRIVATE benchmark_common)

create_single_source_cgal_program("poly_meshing/poly_remesh_feature_preserving.cpp")
target_link_libraries(poly_remesh_feature_preserving PRIVATE benchmark_common)

