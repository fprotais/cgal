#include "statistics.h"
#include "parameters.h"
#include "cgal_common.h"

#include <assert.h>
#include <cmath>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/Polygon_mesh_processing/IO/polygon_mesh_io.h>
#include <CGAL/Polygon_mesh_processing/bbox.h>
#include <CGAL/boost/graph/copy_face_graph.h>
 
#include <CGAL/Mesh_triangulation_3.h>
#include <CGAL/Mesh_complex_3_in_triangulation_3.h>
#include <CGAL/Mesh_criteria_3.h>
 
#include <CGAL/Polyhedron_3.h>
#include <CGAL/Polyhedral_mesh_domain_3.h>
#include <CGAL/make_mesh_3.h>


#include <CGAL/Mesh_smoothing_3/boundary_aware_mesh_smoothing.h>
#include <CGAL/Mesh_smoothing_3/Projectors.h>

using K = CGAL::Exact_predicates_inexact_constructions_kernel;
using Polyhedron = CGAL::Polyhedron_3<K>;
using Mesh_domain = CGAL::Polyhedral_mesh_domain_3<Polyhedron, K>;

// Triangulation
using Tr = CGAL::Mesh_triangulation_3<Mesh_domain>::type;

using C3t3 = CGAL::Mesh_complex_3_in_triangulation_3<Tr>;

// Criteria
using Mesh_criteria = CGAL::Mesh_criteria_3<Tr>;


using Polyhedron_projector = CGAL::Mesh_smoothing_3::Polyhedral_mesh_domain_projector<Mesh_domain>;

namespace params = CGAL::parameters;

int main(int argc, char** argv)
{
    std::string run_name = "poly_remesh";
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
    const CGAL::Bbox_3 bbox = CGAL::Polygon_mesh_processing::bbox(surface);
    const double bbox_diagonal = std::hypot(bbox.xmax() - bbox.xmin(),
                                          bbox.ymax() - bbox.ymin(),
                                          bbox.zmax() - bbox.zmin());
    const double target_size = parameters.bbox_target_sizing * bbox_diagonal;
    Polyhedron polyhedron;
    CGAL::copy_face_graph(surface, polyhedron);

    Mesh_domain domain(polyhedron);
    Mesh_criteria criteria(params::edge_size(0.5 * target_size).
                                 facet_angle(25).facet_size(target_size).facet_distance(0.1 * target_size).
                                 cell_radius_edge_ratio(3).cell_size(target_size));
    
    Timer timer;
    C3t3 c3t3 = CGAL::make_mesh_3<C3t3>(domain, criteria);
    const std::size_t removed_cells = cleanup_non_manifold_cells(c3t3);
    if (parameters.verbose)
        std::cout << "Manifold cleanup removed " << removed_cells << " tetrahedra\n";
    result.init_metrics.computation_time = timer.pop();

    result.nb_faces = c3t3.number_of_facets();
    result.nb_vertices = c3t3.triangulation().number_of_vertices();
    result.nb_tetrahedra = c3t3.number_of_cells();
    
    write_medit(c3t3, prefix + "_initial.mesh");
    write_surface(c3t3, prefix + "_initial.obj");

    result.init_metrics.distance_metrics = compute_distance_metrics(input_name, prefix + "_initial.obj");
    result.init_metrics.gmsh_metrics = compute_gmsh_metrics(prefix + "_initial.mesh");

    // smoothing


    Polyhedron_projector polyhedron_projector(domain);

    const auto status = CGAL::boundary_aware_mesh_smoothing(
        c3t3, polyhedron_projector,
        CGAL::parameters::verbose(parameters.verbose));

    write_medit(c3t3, prefix + "_smoothed.mesh");
    write_surface(c3t3, prefix + "_smoothed.obj");

    result.smoothed_metrics.computation_time = status.total_time;
    result.smoothed_metrics.nb_inverted_elements = status.nb_invalid_elements;
    result.init_metrics.nb_inverted_elements = status.nb_initial_invalid_elements;

    result.smoothed_metrics.distance_metrics = compute_distance_metrics(input_name, prefix + "_smoothed.obj");
    result.smoothed_metrics.gmsh_metrics = compute_gmsh_metrics(prefix + "_smoothed.mesh");

    if (parameters.verbose)
        std::cout << result << std::endl;
    std::cout << result.csv_metric() << std::endl;
    return EXIT_SUCCESS;
}
