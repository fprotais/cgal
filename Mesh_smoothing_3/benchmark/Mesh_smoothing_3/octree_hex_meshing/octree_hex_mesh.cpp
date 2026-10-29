#include "statistics.h"
#include "parameters.h"
#include "cgal_common.h"

#include <assert.h>
#include <cmath>
#include <stdexcept>
#include <vector>

#include <CGAL/Linear_cell_complex_for_combinatorial_map.h>
#include <CGAL/Linear_cell_complex/IO/VTK.h>
#include <CGAL/Linear_cell_complex/IO/lcc_save_load_surfaces.h>
#include <CGAL/generate_hexahedral_mesh_using_two_refinement.h>
#include <CGAL/Hexmeshing/LCC_items.h>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/Polyhedral_mesh_domain_3.h>
#include <CGAL/Polygon_mesh_processing/triangulate_faces.h>
#include <CGAL/Polygon_mesh_processing/orientation.h>
#include <CGAL/Polygon_mesh_processing/IO/polygon_mesh_io.h>
#include <CGAL/Polygon_mesh_processing/bbox.h>

#include <CGAL/Mesh_smoothing_3/mesh_representations.h>
#include <CGAL/Mesh_smoothing_3/Mesh_smoothing_3.h>
#include <CGAL/Mesh_smoothing_3/Projectors.h>


// Preserve the volume fractions produced and corrected by the hexmesher.
using LCC = CGAL::internal::Hexmeshing::LCC;

using K = CGAL::Exact_predicates_inexact_constructions_kernel;
using Point_3 = K::Point_3;

using Surface_mesh = CGAL::Surface_mesh<Point_3>;
using Mesh_domain = CGAL::Polyhedral_mesh_domain_3<Surface_mesh, K>;

namespace params = CGAL::parameters;

// Matches the hexmesher's inner_mark, including its non-manifold corrections.
// Keep this classification after smoothing instead of reclassifying moved cells.
LCC extract_interior(const LCC& lcc) {
    LCC interior(lcc);
    std::vector<LCC::Dart_descriptor> exterior_cells;
    auto cells = interior.one_dart_per_cell<3>();
    for (auto cell = cells.begin(); cell != cells.end(); ++cell) {
        const auto attribute = interior.attribute<3>(cell);
        if (attribute == nullptr)
            throw std::runtime_error("Missing hexmeshing volume-fraction attribute");
        if (!(attribute->info().fraction > 0.5)) exterior_cells.push_back(cell);
    }
    for (auto cell : exterior_cells) interior.remove_cell<3>(cell);
    return interior;
}

void save_surface(LCC const &lcc, std::string const &filename) {
    const LCC interior = extract_interior(lcc);
    CGAL::Surface_mesh<LCC::Point> boundary;
    CGAL::IO::build_manifold_surface_from_lcc(interior, boundary, false);
    // Distance metrics require triangle surfaces, including for hexahedral meshes.
    CGAL::Polygon_mesh_processing::triangulate_faces(boundary);
    CGAL::Polygon_mesh_processing::orient_to_bound_a_volume(boundary);
    CGAL::IO::write_OBJ(filename, boundary, CGAL::parameters::stream_precision(17));
}

void save_interior(const LCC& lcc, const std::string& prefix) {
    const LCC interior = extract_interior(lcc);
    CGAL::IO::write_VTK((prefix + ".vtk").c_str(), interior);
    vtk2medit(prefix + ".vtk", prefix + ".mesh");
}

int main(int argc, char** argv)
{
    std::string run_name = "octree_hex_mesh";
    if (argc != 2) {
        std::cerr << "Usage: " << argv[0]
                  << " input_surface\n";
        std::cout << "Failure" << std::endl;
        return EXIT_FAILURE;
    }

    std::string input_name = argv[1];
    assert(input_name.find_last_of('.') != std::string::npos);
    std::string surface_name =  std::filesystem::path(input_name).stem();

    std::string prefix = run_name + "_" + surface_name;

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


    

    {
        Timer timer;
        LCC lcc=CGAL::generate_hexahedral_mesh_using_two_refinement(
            surface, 
            parameters.octree_cube_cells_per_dim, 
            parameters.octree_nb_levels, 
            CGAL::parameters::use_trimming(false).use_smoothing(true));
        result.init_metrics.computation_time = timer.pop();

        CGAL::IO::write_VTK(std::string(prefix + "_initial.vtk").c_str(), lcc);
        vtk2medit(prefix+"_initial.vtk", prefix + "_initial.mesh");
        save_surface(lcc, prefix + "_initial.obj");
        save_interior(lcc, prefix + "_initial_interior");

        result.init_metrics.distance_metrics = compute_distance_metrics(input_name, prefix + "_initial.obj");
        result.init_metrics.gmsh_metrics = compute_gmsh_metrics(prefix + "_initial.mesh");

    }

    Timer timer;
        LCC lcc=CGAL::generate_hexahedral_mesh_using_two_refinement(
            surface, 
            parameters.octree_cube_cells_per_dim, 
            parameters.octree_nb_levels, 
            CGAL::parameters::use_trimming(false).use_smoothing(false));
    save_interior(lcc, prefix + "_grid_interior");

    result.time_preprocessing = timer.pop();

    // smoothing

    Mesh_domain mesh_domain(surface);
    CGAL::Mesh_smoothing_3::Polyhedral_mesh_domain_projector<Mesh_domain> projector(mesh_domain);

    using LCC_mixed_mesh = CGAL::Mesh_smoothing_3::cgal_types::LCC_mixed_mesh<LCC>;
    using Mixed_mesh = CGAL::Mesh_smoothing_3::helper_structures::Mixed_mesh_wrapper<LCC_mixed_mesh>;

    LCC_mixed_mesh lcc_mixed_mesh(lcc);
    Mixed_mesh mixed_mesh(lcc_mixed_mesh);

    using Lcc_surface = CGAL::Mesh_smoothing_3::cgal_types::LCC_surface_mesh<LCC>;
    Lcc_surface mesh_surface(lcc);

    using Curve_network = CGAL::Mesh_smoothing_3::default_structures::Empty_edge_network<Mixed_mesh::Vertex_descriptor>;

    CGAL::Mesh_smoothing_3::Mesh_smoother<Mixed_mesh, Lcc_surface, Curve_network>
    smoother(mixed_mesh, mesh_surface);

    smoother.set_verbose(parameters.verbose);

    smoother.set_boundary_query([&](std::vector<Point_3> const& pts, Lcc_surface::Surface_patch_index patch_face) {
        auto proj = projector.patch_face_projection_plane(patch_face, pts);
        return std::make_tuple(proj.origin(), proj.vector(), 1.);
    });

    auto status = smoother.run();

    CGAL::IO::write_VTK(std::string(prefix + "_smoothed.vtk").c_str(), lcc);
    vtk2medit(prefix+"_smoothed.vtk", prefix + "_smoothed.mesh");
    save_surface(lcc, prefix + "_smoothed.obj");
    save_interior(lcc, prefix + "_smoothed_interior");

    result.smoothed_metrics.computation_time = status.total_time;
    result.smoothed_metrics.nb_inverted_elements = status.nb_invalid_elements;
    result.init_metrics.nb_inverted_elements = status.nb_initial_invalid_elements;

    result.smoothed_metrics.distance_metrics = compute_distance_metrics(input_name, run_name + "_" + surface_name + "_smoothed.obj");
    result.smoothed_metrics.gmsh_metrics = compute_gmsh_metrics(run_name + "_" + surface_name + "_smoothed.mesh");

    if (parameters.verbose)
        std::cout << result << std::endl;

    std::cout << result.csv_metric() << std::endl;
    return EXIT_SUCCESS;
}
