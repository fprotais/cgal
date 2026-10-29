#include "statistics.h"
#include "parameters.h"
#include "cgal_common.h"

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/alpha_wrap_3.h>
#include <CGAL/Alpha_wrap_3/internal/Triangle_mesh_oracle.h>
#include <CGAL/Delaunay_triangulation_3.h>
#include <CGAL/Robust_weighted_circumcenter_filtered_traits_3.h>
#include <CGAL/Simplicial_mesh_cell_base_3.h>
#include <CGAL/Simplicial_mesh_vertex_base_3.h>
#include <CGAL/Tetrahedral_remeshing/Remeshing_cell_base_3.h>
#include <CGAL/Tetrahedral_remeshing/Remeshing_vertex_base_3.h>
#include <CGAL/tetrahedral_remeshing.h>
#include <CGAL/Mesh_complex_3_in_triangulation_3.h>
#include <CGAL/AABB_tree.h>
#include <CGAL/AABB_traits_3.h>
#include <CGAL/AABB_face_graph_triangle_primitive.h>
#include <CGAL/Polygon_mesh_processing/IO/polygon_mesh_io.h>
#include <CGAL/Polygon_mesh_processing/bbox.h>
#include <CGAL/Polygon_mesh_processing/compute_normal.h>
#include <CGAL/Polygon_mesh_processing/shape_predicates.h>
#include <CGAL/Mesh_smoothing_3/boundary_aware_mesh_smoothing.h>
#include <CGAL/Mesh_smoothing_3/Projectors.h>
#include <CGAL/property_map.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <set>
#include <stdexcept>
#include <string>
#include <vector>

namespace PMP = CGAL::Polygon_mesh_processing;
namespace AW3i = CGAL::Alpha_wraps_3::internal;
namespace params = CGAL::parameters;

using K = CGAL::Exact_predicates_inexact_constructions_kernel;
using Point = K::Point_3;
using Surface_mesh = CGAL::Surface_mesh<Point>;

// The alpha-wrapper's triangulation also needs the bases required by remeshing.
// This is the type stack used in CGAL's volumetric_wrap example.
using Gt = CGAL::Robust_circumcenter_filtered_traits_3<AW3i::Alpha_wrap_AABB_geom_traits<K>>;
using AW_vertex = AW3i::Alpha_wrap_triangulation_vertex_base_3<K>;
using Vertex = CGAL::Simplicial_mesh_vertex_base_3<K, int, int, int, int, AW_vertex>;
using Vb = CGAL::Tetrahedral_remeshing::Remeshing_vertex_base_3<K, Vertex>;
using AW_cell = AW3i::Alpha_wrap_triangulation_cell_base_3<K>;
using Cell = CGAL::Simplicial_mesh_cell_base_3<K, int, int, AW_cell>;
using Cb = CGAL::Tetrahedral_remeshing::Remeshing_cell_base_3<K, Cell>;
using Tds = CGAL::Triangulation_data_structure_3<Vb, Cb>;
using Delaunay = CGAL::Delaunay_triangulation_3<Gt, Tds, CGAL::Fast_location>;
using Tr = CGAL::Triangulation_3<Delaunay::Geom_traits, Delaunay::Triangulation_data_structure>;
using C3t3 = CGAL::Mesh_complex_3_in_triangulation_3<Tr>;
using Primitive = CGAL::AABB_face_graph_triangle_primitive<Surface_mesh>;
using Tree = CGAL::AABB_tree<CGAL::AABB_traits_3<K, Primitive>>;

constexpr int inside = 1;
constexpr int outside = 2;
constexpr int wrap_patch = 1;
constexpr int box_patch = 2;

struct Coarse_sizing {
    const Tree& tree;
    double near_size;
    double far_size;
    double offset;

    template<typename Index>
    double operator()(const Point& p, int, const Index&) const {
        const double distance = std::sqrt(CGAL::to_double(tree.squared_distance(p)));
        // Keep a transition layer near the wrap, then coarsen both volume regions.
        return (std::min)(far_size, near_size + 0.5 * (std::max)(0., distance - offset));
    }
};

class Offset_projector {
public:
    using Tangent_space = CGAL::Mesh_smoothing_3::Tangent_space<K>;

    Offset_projector(const Surface_mesh& surface, const Tree& tree, double offset, double weight)
        : surface(surface), tree(tree), offset(offset), weight(weight) {}

    template<typename Patch_face>
    Tangent_space patch_face_projection_plane(const Patch_face& patch_face,
                                              const std::vector<Point>& points) const {
        if (patch_face.first != wrap_patch)
            return {Point(), K::Vector_3(), CGAL::Mesh_smoothing_3::Projection_weight_mode::NONE};
        const Point center = CGAL::centroid(points.begin(), points.end());
        const auto closest = tree.closest_point_and_primitive(center);
        K::Vector_3 normal = center - closest.first;
        if (normal.squared_length() == 0)
            normal = PMP::compute_face_normal(closest.second, surface);
        normal /= std::sqrt(CGAL::to_double(normal.squared_length()));
        return {closest.first + offset * normal, normal,
                CGAL::Mesh_smoothing_3::Projection_weight_mode::CUSTOM, weight};
    }

    template<typename Curve_edge>
    Tangent_space curve_edge_projection_line(const Curve_edge&,
                                             const std::array<Point, 2>&) const {
        return {Point(), K::Vector_3(), CGAL::Mesh_smoothing_3::Projection_weight_mode::NONE};
    }

private:
    const Surface_mesh& surface;
    const Tree& tree;
    double offset;
    double weight;
};

// Volume outputs contain both regions; surface outputs contain only the wrap.
void write_result(const C3t3& mesh, const std::string& prefix) {
    write_medit(mesh, prefix + ".mesh");
    C3t3 interior(mesh);
    for (auto c : interior.triangulation().finite_cell_handles())
        if (interior.subdomain_index(c) == outside) interior.remove_from_complex(c);
    for (auto f : interior.triangulation().finite_facets())
        if (interior.surface_patch_index(f) == box_patch) interior.remove_from_complex(f);
    write_medit(interior, prefix + "_inside.mesh");
    write_surface(interior, prefix + ".obj");
}

int main(int argc, char** argv) {
    std::string run_name = "offset_meshing";
    if (argc != 2) {
        std::cerr << "Usage: " << argv[0]
                  << " input_surface\n";
        std::cout << "Failure" << std::endl;
        return EXIT_FAILURE;
    }

    std::string input_name = argv[1];
    assert(input_name.find_last_of('.') != std::string::npos);
    std::string surface_name =  std::filesystem::path(input_name).stem();
    std::cout << surface_name << std::endl;

    if (!check_input_surface(input_name)) {
        std::cout << "Failure" << std::endl;
        return EXIT_FAILURE;
    }

    Result_and_statistics result;
    result.run_type = run_name;
    result.case_name = surface_name;
    result.time_preprocessing = 0;

    Input_parameters parameters;
    std::ifstream param_file("parameters.in");
    if (param_file && !(param_file >> parameters)) {
        std::cerr << "Invalid parameters.in\n";
        return EXIT_FAILURE;
    }


    // initial meshing
    CGAL::Surface_mesh<K::Point_3> surface;
    CGAL::Polygon_mesh_processing::IO::read_polygon_mesh(input_name, surface);
    const CGAL::Bbox_3 bbox = CGAL::Polygon_mesh_processing::bbox(surface);
    const double bbox_diagonal = std::hypot(bbox.xmax() - bbox.xmin(),
                                          bbox.ymax() - bbox.ymin(),
                                          bbox.zmax() - bbox.zmin());
    const double target_size = parameters.bbox_target_sizing * bbox_diagonal;
    const double delta = parameters.offset * bbox_diagonal;

    // edge length of regular tetrahedron with circumradius alpha
    const double alpha = target_size / 1.6329931618554521; // sqrt(8/3)

    Timer preprocessing;
    Tree tree(faces(surface).first, faces(surface).second, surface);
    tree.accelerate_distance_queries();
    result.time_preprocessing = preprocessing.pop();


    Timer meshing;
    AW3i::Triangle_mesh_oracle<K> oracle(target_size);
    oracle.add_triangle_mesh(surface);
    AW3i::Alpha_wrapper_3<decltype(oracle), Delaunay> wrapper(oracle);
    Surface_mesh wrap;
    wrapper(alpha, delta, wrap);
    Tr tr(static_cast<const Tr&>(wrapper.triangulation()));

    // Label every finite cell, so both the inside and outside are remeshed.
    for (auto c : tr.all_cell_handles()) {
        c->set_subdomain_index(tr.is_infinite(c) ? 0 : (c->is_outside() ? outside : inside));
        for (int i = 0; i < 4; ++i) c->set_surface_patch_index(i, 0);
    }
    for (auto v : tr.finite_vertex_handles()) v->set_dimension(3);
    std::set<Tr::Facet> locked_faces;
    std::set<Tr::Vertex_handle> locked_vertices;
    for (auto f : tr.finite_facets()) {
        const auto mirror = tr.mirror_facet(f);
        const int a = f.first->subdomain_index();
        const int b = mirror.first->subdomain_index();
        if (a == b) continue;
        const int patch = (a == 0 || b == 0) ? box_patch : wrap_patch;
        f.first->set_surface_patch_index(f.second, patch);
        mirror.first->set_surface_patch_index(mirror.second, patch);
        for (auto v : tr.vertices(f)) v->set_dimension(2);
        if (patch == wrap_patch) {
            locked_faces.insert(f);
            locked_faces.insert(mirror);
            for (auto v : tr.vertices(f)) locked_vertices.insert(v);
        }
    }

    const Coarse_sizing sizing{tree, std::sqrt(8. / 3.) * target_size,
                                (std::max)(4. * target_size, 0.15 * bbox_diagonal), delta};
    // Constrain only the wrap; the external box may be remeshed.
    CGAL::tetrahedral_isotropic_remeshing(tr, sizing,
        params::remesh_boundaries(true)
               .facet_is_constrained_map(CGAL::Boolean_property_map<std::set<Tr::Facet>>(locked_faces))
               .vertex_is_constrained_map(CGAL::Boolean_property_map<std::set<Tr::Vertex_handle>>(locked_vertices))
               .number_of_iterations(3));
    // Remeshing marks constrained vertices as corners. Release these temporary
    // locks before constructing the complex used for offset smoothing.
    for (auto v : locked_vertices) {
        v->set_dimension(2);
        v->set_index(Tr::Vertex::Index(wrap_patch));
    }
    C3t3 c3t3;
    c3t3.set_triangulation(std::move(tr));
    result.init_metrics.computation_time = meshing.pop();
    result.nb_vertices = c3t3.triangulation().number_of_vertices();
    result.nb_tetrahedra = c3t3.number_of_cells();
    result.nb_faces = 0;
    for (auto f : c3t3.facets_in_complex())
        if (c3t3.surface_patch_index(f) == wrap_patch) ++result.nb_faces;

    const std::string prefix = result.run_type + "_" + surface_name;
    write_result(c3t3, prefix + "_initial");
    result.init_metrics.distance_metrics = compute_offset_distance_metrics(input_name, prefix + "_initial.obj", delta);
    result.init_metrics.gmsh_metrics = compute_gmsh_metrics(prefix + "_initial.mesh");
    
    // smoothing
    // Release the wrap for offset projection, but anchor the external box.
    std::set<Tr::Vertex_handle> box_vertices;
    for (auto f : c3t3.facets_in_complex())
        if (c3t3.surface_patch_index(f) == box_patch)
            for (auto v : c3t3.triangulation().vertices(f)) box_vertices.insert(v);
    const Offset_projector projector(surface, tree, delta, parameters.boundary_weight);
    const auto status = CGAL::boundary_aware_mesh_smoothing(c3t3, projector,
        params::vertex_is_constrained_map(CGAL::Boolean_property_map<std::set<Tr::Vertex_handle>>(box_vertices))
                .verbose(parameters.verbose));
                
    result.smoothed_metrics.computation_time = status.total_time;
    result.smoothed_metrics.nb_inverted_elements = status.nb_invalid_elements;
    result.init_metrics.nb_inverted_elements = status.nb_initial_invalid_elements;

    write_result(c3t3, prefix + "_smoothed");
    result.smoothed_metrics.distance_metrics = compute_offset_distance_metrics(input_name, prefix + "_smoothed.obj", delta);
    result.smoothed_metrics.gmsh_metrics = compute_gmsh_metrics(prefix + "_smoothed.mesh");

    if (parameters.verbose)
        std::cout << result << std::endl;
    std::cout << result.csv_metric() << std::endl;
    return EXIT_SUCCESS;
}
