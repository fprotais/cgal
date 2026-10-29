#include "statistics.h"
#include "parameters.h"
#include "cgal_common.h"

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/Polygon_mesh_processing/IO/polygon_mesh_io.h>
#include <CGAL/Polygon_mesh_processing/bbox.h>
#include <CGAL/AABB_tree.h>
#include <CGAL/AABB_traits_3.h>
#include <CGAL/AABB_face_graph_triangle_primitive.h>
#include <CGAL/Side_of_triangle_mesh.h>
#include <CGAL/Labeled_mesh_domain_3.h>
#include <CGAL/Mesh_triangulation_3.h>
#include <CGAL/Mesh_complex_3_in_triangulation_3.h>
#include <CGAL/Mesh_criteria_3.h>
#include <CGAL/make_mesh_3.h>
#include <CGAL/squared_distance_3.h>
#include <CGAL/Mesh_smoothing_3/boundary_aware_mesh_smoothing.h>
#include <CGAL/Mesh_smoothing_3/Projectors.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <limits>
#include <optional>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

using Kernel = CGAL::Exact_predicates_inexact_constructions_kernel;
using Point = Kernel::Point_3;
using Vector = Kernel::Vector_3;
using Surface_mesh = CGAL::Surface_mesh<Point>;
namespace PMP = CGAL::Polygon_mesh_processing;
namespace params = CGAL::parameters;

using Primitive = CGAL::AABB_face_graph_triangle_primitive<Surface_mesh>;
using AABB_traits = CGAL::AABB_traits_3<Kernel, Primitive>;
using AABB_tree = CGAL::AABB_tree<AABB_traits>;
using Side = CGAL::Side_of_triangle_mesh<Surface_mesh, Kernel>;
using Domain = CGAL::Labeled_mesh_domain_3<Kernel>;
using Triangulation = CGAL::Mesh_triangulation_3<Domain>::type;
using C3t3 = CGAL::Mesh_complex_3_in_triangulation_3<Triangulation>;
using Criteria = CGAL::Mesh_criteria_3<Triangulation>;

class SDF_grid {
public:
    struct Sample {
        double value;
        Vector gradient;
    };

    SDF_grid(const Surface_mesh& surface, std::size_t n)
        : n_(n)
    {
        if (n_ < 4 || n_ > 512)
            throw std::runtime_error("SDF grid_size must be between 4 and 512");
        data_.resize(n_ * n_ * n_);
        std::array<double, 3> lo = {
            std::numeric_limits<double>::max(),
            std::numeric_limits<double>::max(),
            std::numeric_limits<double>::max()
        };
        std::array<double, 3> hi = {
            std::numeric_limits<double>::lowest(),
            std::numeric_limits<double>::lowest(),
            std::numeric_limits<double>::lowest()
        };
        for (auto v : surface.vertices()) {
            const Point& p = surface.point(v);
            const std::array<double, 3> x = {p.x(), p.y(), p.z()};
            for (int a = 0; a < 3; ++a) {
                lo[a] = std::min(lo[a], x[a]);
                hi[a] = std::max(hi[a], x[a]);
            }
        }

        const double span = std::max({hi[0] - lo[0], hi[1] - lo[1], hi[2] - lo[2]});

        // A cubic grid, 15% padding on the largest bounding-box dimension.
        const double length = 1.3 * span;
        spacing_ = length / double(n_ - 1);
        for (int a = 0; a < 3; ++a) {
            minimum_[a] = 0.5 * (lo[a] + hi[a] - length);
            maximum_[a] = minimum_[a] + length;
        }

        // Both distance and bounded-side queries share this AABB tree.
        const auto face_range = faces(surface);
        AABB_tree tree(face_range.begin(), face_range.end(), surface);
        tree.accelerate_distance_queries();
        const Side side(tree);

        double most_negative = 0.0;
        for (std::size_t i = 0; i < n_; ++i) {
            for (std::size_t j = 0; j < n_; ++j) {
                for (std::size_t k = 0; k < n_; ++k) {
                    const Point p = grid_point(i, j, k);
                    double distance = std::sqrt(CGAL::to_double(tree.squared_distance(p)));
                    const CGAL::Bounded_side location = side(p);
                    if (location == CGAL::ON_BOUNDED_SIDE) distance = -distance;
                    if (location == CGAL::ON_BOUNDARY) distance = 0.0;
                    data_[index(i, j, k)] = distance;
                    if (distance < most_negative) {
                        most_negative = distance;
                        interior_ = p;
                    }
                }
            }
        }
    }

    double spacing() const { return spacing_; }
    Point interior_point() const { return *interior_; }

    // Sphere centered inside the sampled negative region, enclosing the cube.
    Kernel::Sphere_3 bounding_sphere() const
    {
        double r2 = 0.0;
        for (int i = 0; i < 2; ++i)
            for (int j = 0; j < 2; ++j)
                for (int k = 0; k < 2; ++k) {
                    const Point corner(i ? maximum_[0] : minimum_[0],
                                       j ? maximum_[1] : minimum_[1],
                                       k ? maximum_[2] : minimum_[2]);
                    r2 = std::max(r2, CGAL::to_double(
                                              CGAL::squared_distance(*interior_, corner)));
                }
        return Kernel::Sphere_3(*interior_, 1.01 * r2);
    }

    // Samples a continuous extension of the grid function to all of R^3.
    // Outside the cube: f(p) = f(clamp(p)) + distance(p, clamp(p)).
    Sample evaluate(const Point& p) const
    {
        const std::array<double, 3> x = {p.x(), p.y(), p.z()};
        std::array<double, 3> q, delta;
        double outside_sq = 0.0;
        for (int a = 0; a < 3; ++a) {
            q[a] = std::clamp(x[a], minimum_[a], maximum_[a]);
            delta[a] = x[a] - q[a];
            outside_sq += delta[a] * delta[a];
        }

        Sample result = interpolate(q);
        if (outside_sq > 0.0) {
            const double outside = std::sqrt(outside_sq);
            result.value += outside;
            std::array<double, 3> g = {
                result.gradient.x(), result.gradient.y(), result.gradient.z()
            };
            // In clamped directions, the clamped point does not move with p.
            for (int a = 0; a < 3; ++a)
                if (delta[a] != 0.0) g[a] = delta[a] / outside;
            result.gradient = Vector(g[0], g[1], g[2]);
        }
        return result;
    }

    double value(const Point& p) const { return evaluate(p).value; }
    Vector gradient(const Point& p) const { return evaluate(p).gradient; }

private:
    std::size_t index(std::size_t i, std::size_t j, std::size_t k) const
    { return (i * n_ + j) * n_ + k; }

    Point grid_point(std::size_t i, std::size_t j, std::size_t k) const
    {
        return Point(minimum_[0] + i * spacing_,
                     minimum_[1] + j * spacing_,
                     minimum_[2] + k * spacing_);
    }

    // Trilinear value and its exact (cellwise) partial derivatives.
    Sample interpolate(const std::array<double, 3>& q) const
    {
        std::array<std::size_t, 3> cell;
        std::array<double, 3> t;
        for (int a = 0; a < 3; ++a) {
            const double u = std::clamp((q[a] - minimum_[a]) / spacing_,
                                        0.0, double(n_ - 1));
            cell[a] = std::min(std::size_t(std::floor(u)), n_ - 2);
            t[a] = u - double(cell[a]);
        }
        double f = 0.0;
        std::array<double, 3> grad = {0.0, 0.0, 0.0};
        for (int a = 0; a < 2; ++a)
            for (int b = 0; b < 2; ++b)
                for (int c = 0; c < 2; ++c) {
                    const double v = data_[index(cell[0] + a,
                                                 cell[1] + b,
                                                 cell[2] + c)];
                    const double wx = a ? t[0] : 1.0 - t[0];
                    const double wy = b ? t[1] : 1.0 - t[1];
                    const double wz = c ? t[2] : 1.0 - t[2];
                    f += v * wx * wy * wz;
                    grad[0] += v * (a ? 1.0 : -1.0) * wy * wz / spacing_;
                    grad[1] += v * wx * (b ? 1.0 : -1.0) * wz / spacing_;
                    grad[2] += v * wx * wy * (c ? 1.0 : -1.0) / spacing_;
                }
        return {f, Vector(grad[0], grad[1], grad[2])};
    }

    std::size_t n_;
    std::vector<double> data_;
    double spacing_ = 0.0;
    std::array<double, 3> minimum_{}, maximum_{};
    std::optional<Point> interior_;
};

// Area-weighted mean absolute SDF value and normal angle (degrees).
// Feature and corner distances are zero for this smooth implicit representation.
Distance_metrics compute_sdf_distance_metrics(const C3t3& c3t3, const SDF_grid& sdf) {
    Distance_metrics metrics{};
    Surface_mesh boundary;
    CGAL::facets_in_complex_3_to_triangle_mesh(c3t3, boundary);

    const double rad_to_deg = 180. / std::acos(-1.);

    double total_area = 0.0;
    double distance_error = 0.0;
    double normal_error = 0.0;
    double max_error = 0.0;

    for (auto f : boundary.faces()) {
        auto h0 = boundary.halfedge(f);
        auto h1 = boundary.next(h0);
        auto h2 = boundary.next(h1);

        const Point& a = boundary.point(boundary.source(h0));
        const Point& b = boundary.point(boundary.source(h1));
        const Point& c = boundary.point(boundary.source(h2));

        const Vector n = CGAL::cross_product(b - a, c - a);
        const double n_length =
            std::sqrt(CGAL::to_double(n.squared_length()));
        const double area = 0.5 * n_length;

        total_area += area;

        // Subdivide triangles to sample at roughly the SDF grid spacing.
        const double max_edge = std::sqrt(std::max({
            CGAL::to_double(CGAL::squared_distance(a, b)),
            CGAL::to_double(CGAL::squared_distance(b, c)),
            CGAL::to_double(CGAL::squared_distance(c, a))
        }));

        const int m = std::max(1,
            static_cast<int>(std::ceil(max_edge / sdf.spacing())));
        const double sample_area = area / (double(m) * double(m));

        // Evaluate the centroid of each small triangle.
        auto sample = [&](double u, double v) {
            const Point p = a + (b - a) * u + (c - a) * v;
            const auto s = sdf.evaluate(p);

            distance_error += sample_area * std::abs(s.value);
            max_error = std::max(max_error, std::abs(s.value));

            const double g_length = std::sqrt(
                CGAL::to_double(s.gradient.squared_length()));

            double angle = 90.0; // Undefined gradient: penalty
            if (g_length > 0.0) {
                const double cosine = std::clamp(
                    CGAL::to_double(n * s.gradient) /
                        (n_length * g_length),
                    -1.0, 1.0);
                angle = std::acos(cosine) * rad_to_deg;
            }

            normal_error += sample_area * angle;
        };

        for (int i = 0; i < m; ++i) {
            for (int j = 0; j < m - i; ++j) {
                sample((i + 1.0 / 3.0) / m,
                       (j + 1.0 / 3.0) / m);

                if (i + j + 1 < m)
                    sample((i + 2.0 / 3.0) / m,
                           (j + 2.0 / 3.0) / m);
            }
        }

        // Include vertices and edge midpoints in the maximum.
        for (const Point& p : {
                 a, b, c,
                 a + (b - a) * 0.5,
                 b + (c - b) * 0.5,
                 c + (a - c) * 0.5
             })
            max_error = std::max(max_error, std::abs(sdf.value(p)));
    }

    metrics.l2_distance = distance_error / total_area;
    metrics.normal_deviation = normal_error / total_area;
    // Sampled maximum distance to the SDF zero set, rather than to the input mesh.
    metrics.hausdorff_distance = max_error;
    metrics.features_hausdorff_distance = 0.;
    metrics.features_l2_distance = 0.;
    metrics.corner_hausdorff_distance = 0.;
    return metrics;
}

int main(int argc, char** argv) {
    const std::string run_name = "sdf_meshing";
    if (argc != 2) {
        std::cerr << "Usage: " << argv[0] << " input_surface\n";
        std::cout << "Failure" << std::endl;
        return EXIT_FAILURE;
    }

    const std::string input_name = argv[1];
    const std::string surface_name = std::filesystem::path(input_name).stem().string();
    const std::string prefix = run_name + "_" + surface_name;
    if (!check_input_surface(input_name)) {
        std::cout << "Failure" << std::endl;
        return EXIT_FAILURE;
    }

    Result_and_statistics result{};
    result.run_type = run_name;
    result.case_name = surface_name;

    Input_parameters parameters;
    std::ifstream param_file("parameters.in");
    if (param_file && !(param_file >> parameters)) {
        std::cerr << "Invalid parameters.in\n";
        return EXIT_FAILURE;
    }

    Surface_mesh surface;
    CGAL::Polygon_mesh_processing::IO::read_polygon_mesh(input_name, surface);
    const auto bbox = PMP::bbox(surface);
    const double bbox_diagonal = std::hypot(bbox.xmax() - bbox.xmin(),
                                            bbox.ymax() - bbox.ymin(),
                                            bbox.zmax() - bbox.zmin());
    const double target_size = parameters.bbox_target_sizing * bbox_diagonal;

    Timer preprocessing;
    const SDF_grid sdf(surface, parameters.grid_size);
    result.time_preprocessing = preprocessing.pop();

    // The sampled SDF defines both the meshing domain and the projection target.
    const auto implicit = [&sdf](const Point& p) { return sdf.value(p); };
    const Domain domain = Domain::create_implicit_mesh_domain(
        implicit, sdf.bounding_sphere(), params::relative_error_bound(1e-5));
    Criteria criteria(params::edge_size(0.5 * target_size).
                                 facet_angle(25).facet_size(target_size).facet_distance(0.1 * target_size).
                                 cell_radius_edge_ratio(3).cell_size(target_size));
    Timer timer;
    C3t3 c3t3 = CGAL::make_mesh_3<C3t3>(domain, criteria);
    result.init_metrics.computation_time = timer.pop();
    result.nb_faces = c3t3.number_of_facets();
    result.nb_vertices = c3t3.triangulation().number_of_vertices();
    result.nb_tetrahedra = c3t3.number_of_cells();

    write_medit(c3t3, prefix + "_initial.mesh");
    write_surface(c3t3, prefix + "_initial.obj");
    result.init_metrics.distance_metrics = compute_sdf_distance_metrics(c3t3, sdf);
    result.init_metrics.gmsh_metrics = compute_gmsh_metrics(prefix + "_initial.mesh");

    // Projection using the sampled SDF and its analytical gradient.
    const auto projection_function = [&sdf](const Point& p) {
        const auto sample = sdf.evaluate(p);
        return std::make_pair(sample.value, sample.gradient);
    };
    using Projector = CGAL::Mesh_smoothing_3::Signed_distance_function_projector<
        Kernel, decltype(projection_function)>;
    const Projector projector(projection_function, 10, 1e-4 * target_size);
    const auto status = CGAL::boundary_aware_mesh_smoothing(
        c3t3, projector, params::verbose(parameters.verbose));

    write_medit(c3t3, prefix + "_smoothed.mesh");
    write_surface(c3t3, prefix + "_smoothed.obj");
    result.smoothed_metrics.computation_time = status.total_time;
    result.smoothed_metrics.nb_inverted_elements = status.nb_invalid_elements;
    result.init_metrics.nb_inverted_elements = status.nb_initial_invalid_elements;
    result.smoothed_metrics.distance_metrics = compute_sdf_distance_metrics(c3t3, sdf);
    result.smoothed_metrics.gmsh_metrics = compute_gmsh_metrics(prefix + "_smoothed.mesh");

    if (parameters.verbose)
        std::cout << result << std::endl;
    std::cout << result.csv_metric() << std::endl;
    return EXIT_SUCCESS;

}
