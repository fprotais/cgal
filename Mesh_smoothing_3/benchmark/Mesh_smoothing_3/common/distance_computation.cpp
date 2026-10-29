#include "distance_computation.h"

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/AABB_tree.h>
#include <CGAL/AABB_traits_3.h>
#include <CGAL/AABB_face_graph_triangle_primitive.h>
#include <CGAL/Polyhedral_mesh_domain_with_features_3.h>
#include <CGAL/Polygon_mesh_processing/IO/polygon_mesh_io.h>
#include <CGAL/Polygon_mesh_processing/bbox.h>
#include <CGAL/Polygon_mesh_processing/compute_normal.h>
#include <CGAL/Polygon_mesh_processing/distance.h>
#include <CGAL/Polygon_mesh_processing/measure.h>
#include <CGAL/Polygon_mesh_processing/shape_predicates.h>
#include <CGAL/boost/graph/copy_face_graph.h>
#include <CGAL/point_generators_3.h>
#include <CGAL/squared_distance_3.h>

#include <algorithm>
#include <cmath>
#include <iterator>
#include <stdexcept>
#include <utility>
#include <vector>

namespace {

namespace PMP = CGAL::Polygon_mesh_processing;
using Kernel = CGAL::Exact_predicates_inexact_constructions_kernel;
using Point = Kernel::Point_3;
using Vector = Kernel::Vector_3;
using Surface_mesh = CGAL::Surface_mesh<Point>;
using Face = Surface_mesh::Face_index;
using Face_normals = Surface_mesh::Property_map<Face, Vector>;
using Primitive = CGAL::AABB_face_graph_triangle_primitive<Surface_mesh>;
using Tree = CGAL::AABB_tree<CGAL::AABB_traits_3<Kernel, Primitive>>;
using Mesh_domain = CGAL::Polyhedral_mesh_domain_with_features_3<Kernel>;

constexpr std::size_t surface_sample_count = 20000;
constexpr std::size_t feature_sample_count = 20000;
constexpr double relative_hausdorff_error = 1e-4;

double bbox_diagonal(const Surface_mesh& mesh) {
    const CGAL::Bbox_3 bbox = PMP::bbox(mesh);
    const double diagonal = std::hypot(bbox.xmax() - bbox.xmin(),
                                      bbox.ymax() - bbox.ymin(),
                                      bbox.zmax() - bbox.zmin());
    return diagonal;
}

// Report every length as a percentage of the reference bbox diagonal.
// Normal deviation is an angle and remains in degrees.
Distance_metrics normalize_distances(Distance_metrics metrics, double diagonal) {
    const double factor = 100. / diagonal;
    metrics.hausdorff_distance *= factor;
    metrics.l2_distance *= factor;
    metrics.features_hausdorff_distance *= factor;
    metrics.features_l2_distance *= factor;
    metrics.corner_hausdorff_distance *= factor;
    return metrics;
}

struct Directed_metrics {
    double max_distance = 0;
    double mean_distance = 0;
    double mean_normal_deviation = 0;
};

Directed_metrics sample_surface(const Surface_mesh& source, const Tree& target,
                               const Face_normals& source_normals,
                               const Face_normals& target_normals,
                               bool compare_normals, double offset = 0.) {
    Directed_metrics result;
    const double total_area = CGAL::to_double(PMP::area(source));
    const double radians_to_degrees = 180. / std::acos(-1.);
    CGAL::Random random(0);

    for (Face f : source.faces()) {
        const double area_fraction = CGAL::to_double(PMP::face_area(f, source)) / total_area;
        const std::size_t count = (std::max)(std::size_t(1),
            static_cast<std::size_t>(std::ceil(surface_sample_count * area_fraction)));
        const double sample_weight = area_fraction / count;
        const auto h = source.halfedge(f);
        CGAL::Random_points_in_triangle_3<Point> samples(
            source.point(source.source(h)), source.point(source.target(h)),
            source.point(source.target(source.next(h))), random);

        for (std::size_t i = 0; i < count; ++i, ++samples) {
            const Point& p = *samples;
            const auto closest = target.closest_point_and_primitive(p);
            const double distance = std::abs(
                std::sqrt(CGAL::to_double(CGAL::squared_distance(p, closest.first))) - offset);
            result.mean_distance += sample_weight * distance;
            result.max_distance = (std::max)(result.max_distance, distance);
            if (compare_normals) {
                const Vector& n = source_normals[f];
                const Vector& m = target_normals[closest.second];
                const double sine = std::sqrt(CGAL::to_double(CGAL::cross_product(n, m).squared_length()));
                const double cosine = CGAL::to_double(n * m);
                result.mean_normal_deviation += sample_weight
                    * std::atan2(sine, cosine) * radians_to_degrees;
            }
        }
    }
    return result;
}

void compute_feature_metrics(const Surface_mesh& reference, const Tree& generated,
                             Distance_metrics& metrics, double offset = 0.) {
    Mesh_domain::Polyhedron polyhedron;
    CGAL::copy_face_graph(reference, polyhedron);
    Mesh_domain domain(polyhedron);
    domain.detect_features();

    // Sample detected edges directly, avoiding repeated walks along polylines.
    std::vector<Kernel::Segment_3> segments;
    double total_length = 0;
    for (const auto& surface : domain.polyhedra()) {
        const auto is_feature = get(CGAL::edge_is_feature, surface);
        for (auto e : edges(surface)) {
            if (!get(is_feature, e)) continue;
            const auto h = halfedge(e, surface);
            const auto points = get(CGAL::vertex_point, surface);
            segments.emplace_back(get(points, source(h, surface)), get(points, target(h, surface)));
            total_length += std::sqrt(CGAL::to_double(segments.back().squared_length()));
        }
    }

    const auto distance = [&generated, offset](const Point& point) {
        return std::abs(std::sqrt(CGAL::to_double(generated.squared_distance(point))) - offset);
    };
    for (const auto& segment : segments) {
        const double length = std::sqrt(CGAL::to_double(segment.squared_length()));
        const std::size_t count = (std::max)(std::size_t(1),
            static_cast<std::size_t>(std::ceil(feature_sample_count * length / total_length)));
        const double weight = length / total_length / count;
        // Endpoints contribute to the maximum, not to the length-weighted mean.
        metrics.features_hausdorff_distance = (std::max)({metrics.features_hausdorff_distance,
            distance(segment.source()), distance(segment.target())});
        for (std::size_t i = 0; i < count; ++i) {
            const double t = (i + 0.5) / count;
            const Point p = segment.source() + t * (segment.target() - segment.source());
            const double d = distance(p);
            metrics.features_l2_distance += weight * d;
            metrics.features_hausdorff_distance = (std::max)(metrics.features_hausdorff_distance, d);
        }
    }

    // Use domain corners: junctions and sharp turns, rather than every degree-2
    // vertex along a tessellated feature curve.
    std::vector<std::pair<Mesh_domain::Corner_index, Point>> corners;
    domain.get_corners(std::back_inserter(corners));
    for (const auto& corner : corners)
        metrics.corner_hausdorff_distance += distance(corner.second) / corners.size();
}

} // namespace

Distance_metrics compute_distance_metrics(const std::string& reference_surface,
                                         const std::string& generated_surface) {
    Surface_mesh reference;
    PMP::IO::read_polygon_mesh(reference_surface, reference);
    Surface_mesh generated;
    PMP::IO::read_polygon_mesh(generated_surface, generated);
    Tree reference_tree(faces(reference).first, faces(reference).second, reference);
    Tree generated_tree(faces(generated).first, faces(generated).second, generated);
    reference_tree.accelerate_distance_queries();
    generated_tree.accelerate_distance_queries();

    const auto reference_normals = reference.add_property_map<Face, Vector>("f:normal", CGAL::NULL_VECTOR).first;
    const auto generated_normals = generated.add_property_map<Face, Vector>("f:normal", CGAL::NULL_VECTOR).first;
    PMP::compute_face_normals(reference, reference_normals);
    PMP::compute_face_normals(generated, generated_normals);

    const auto forward = sample_surface(reference, generated_tree, reference_normals, generated_normals, false);
    const auto backward = sample_surface(generated, reference_tree, generated_normals, reference_normals, true);

    const double diagonal = bbox_diagonal(reference);

    Distance_metrics metrics{};
    metrics.hausdorff_distance = PMP::bounded_error_symmetric_Hausdorff_distance<CGAL::Sequential_tag>(
        reference, generated, relative_hausdorff_error * diagonal);
    metrics.l2_distance = 0.5 * (forward.mean_distance + backward.mean_distance);
    metrics.normal_deviation = backward.mean_normal_deviation;
    compute_feature_metrics(reference, generated_tree, metrics);
    return normalize_distances(metrics, diagonal);
}

Distance_metrics compute_offset_distance_metrics(const std::string& reference_surface,
                                                const std::string& generated_surface,
                                                double offset) {
    if (offset == 0.)
        return compute_distance_metrics(reference_surface, generated_surface);

    Surface_mesh reference, generated;
    PMP::IO::read_polygon_mesh(reference_surface, reference);
    PMP::IO::read_polygon_mesh(generated_surface, generated);
    const double diagonal = bbox_diagonal(reference);

    Tree reference_tree(faces(reference).first, faces(reference).second, reference);
    Tree generated_tree(faces(generated).first, faces(generated).second, generated);
    reference_tree.accelerate_distance_queries();
    generated_tree.accelerate_distance_queries();
    const auto reference_normals = reference.add_property_map<Face, Vector>("f:normal", CGAL::NULL_VECTOR).first;
    const auto generated_normals = generated.add_property_map<Face, Vector>("f:normal", CGAL::NULL_VECTOR).first;
    PMP::compute_face_normals(reference, reference_normals);
    PMP::compute_face_normals(generated, generated_normals);

    // Same query as Offset_projector::patch_face_projection_plane in offset_meshing.
    const auto project = [&](const Point& p) {
        const auto closest = reference_tree.closest_point_and_primitive(p);
        Vector normal = p - closest.first;
        if (normal.squared_length() == 0)
            normal = reference_normals[closest.second];
        normal /= std::sqrt(CGAL::to_double(normal.squared_length()));
        return std::make_pair(closest.first + offset * normal, normal);
    };

    const auto forward = sample_surface(reference, generated_tree, reference_normals,
                                        generated_normals, false, offset);
    Directed_metrics backward;
    const double total_area = CGAL::to_double(PMP::area(generated));
    const double radians_to_degrees = 180. / std::acos(-1.);
    CGAL::Random random(0);
    for (Face f : generated.faces()) {
        const double area_fraction = CGAL::to_double(PMP::face_area(f, generated)) / total_area;
        const std::size_t count = (std::max)(std::size_t(1),
            static_cast<std::size_t>(std::ceil(surface_sample_count * area_fraction)));
        const double weight = area_fraction / count;
        const auto h = generated.halfedge(f);
        CGAL::Random_points_in_triangle_3<Point> samples(
            generated.point(generated.source(h)), generated.point(generated.target(h)),
            generated.point(generated.target(generated.next(h))), random);
        for (std::size_t i = 0; i < count; ++i, ++samples) {
            const Point p = *samples;
            const auto projection = project(p);
            const Point& projected = projection.first;
            const double distance = std::sqrt(CGAL::to_double(CGAL::squared_distance(p, projected)));
            backward.mean_distance += weight * distance;
            backward.max_distance = (std::max)(backward.max_distance, distance);
            const Vector& normal = projection.second;
            const Vector n = generated_normals[f];
            const double sine = std::sqrt(CGAL::to_double(CGAL::cross_product(n, normal).squared_length()));
            backward.mean_normal_deviation += weight
                * std::atan2(sine, CGAL::to_double(n * normal)) * radians_to_degrees;
        }
    }

    Distance_metrics metrics{};
    // The offset is implicit: use a sampled symmetric maximum, supplemented
    // with all vertices, rather than PMP's explicit-mesh Hausdorff routine.
    metrics.hausdorff_distance = (std::max)(forward.max_distance, backward.max_distance);
    for (auto v : reference.vertices())
        metrics.hausdorff_distance = (std::max)(metrics.hausdorff_distance,
            std::abs(std::sqrt(CGAL::to_double(generated_tree.squared_distance(reference.point(v)))) - offset));
    for (auto v : generated.vertices()) {
        const Point& p = generated.point(v);
        metrics.hausdorff_distance = (std::max)(metrics.hausdorff_distance,
            std::sqrt(CGAL::to_double(CGAL::squared_distance(p, project(p).first))));
    }

    metrics.l2_distance = 0.5 * (forward.mean_distance + backward.mean_distance);
    metrics.normal_deviation = backward.mean_normal_deviation;
    compute_feature_metrics(reference, generated_tree, metrics, offset);
    return normalize_distances(metrics, diagonal);
}


