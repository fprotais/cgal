#ifndef CGAL_COMMON_H
#define CGAL_COMMON_H

#include <string>
#include <fstream>
#include <algorithm>
#include <cmath>
#include <iterator>
#include <set>
#include <stdexcept>
#include <vector>

#include <CGAL/IO/File_medit.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/boost/graph/IO/OBJ.h>
#include <CGAL/facets_in_complex_3_to_triangle_mesh.h>

bool check_input_surface(std::string input);

// Erode a single-region, single-patch complex until its volume boundary is
// manifold. Features and multiple material/patch indices are deliberately not
// supported: removing cells would require reconstructing those constraints.
// Returns the number of removed tetrahedra. No triangulation cells are deleted.
template <typename C3t3>
std::size_t cleanup_non_manifold_cells(C3t3& c3t3) {
    using Cell = typename C3t3::Cell_handle;
    using Vertex = typename C3t3::Vertex_handle;
    using Patch = typename C3t3::Surface_patch_index;
    using Subdomain = typename C3t3::Subdomain_index;
    auto& tr = c3t3.triangulation();
    if (tr.dimension() != 3 || c3t3.number_of_cells() == 0)
        throw std::invalid_argument("Cleanup requires a nonempty 3D volume complex");
    if (c3t3.number_of_edges_in_complex() || c3t3.number_of_vertices_in_complex())
        throw std::invalid_argument("Cleanup does not support protected features");
    const Subdomain subdomain = c3t3.subdomain_index(c3t3.cells_in_complex_begin());
    for (Cell c : c3t3.cells_in_complex())
        if (c3t3.subdomain_index(c) != subdomain)
            throw std::invalid_argument("Cleanup requires a single volume region");
    Patch patch{};
    for (auto f : c3t3.facets_in_complex()) {
        const Patch current = c3t3.surface_patch_index(f);
        if (patch == Patch()) patch = current;
        else if (current != patch)
            throw std::invalid_argument("Cleanup requires a single surface patch");
    }
    if (patch == Patch())
        throw std::invalid_argument("Cleanup requires a surface patch index");

    const auto inside = [&](Cell c) {
        return !tr.is_infinite(c) && c3t3.is_in_complex(c);
    };
    // As in Alpha_wrapper_3::is_non_manifold(), flood both labels in the
    // vertex umbrella. Include infinite cells to detect pinching at the hull.
    const auto is_non_manifold = [&](Vertex v) {
        std::vector<Cell> incident;
        tr.incident_cells(v, std::back_inserter(incident));
        std::set<Cell> visited;
        unsigned inside_components = 0, outside_components = 0;
        for (Cell seed : incident) {
            if (!visited.insert(seed).second) continue;
            const bool label = inside(seed);
            if (++(label ? inside_components : outside_components) > 1)
                return true;
            std::vector<Cell> pending{seed};
            while (!pending.empty()) {
                Cell c = pending.back();
                pending.pop_back();
                const int vertex_index = c->index(v);
                for (int i = 0; i < 4; ++i) {
                    if (i == vertex_index) continue;
                    Cell neighbor = c->neighbor(i);
                    if (inside(neighbor) == label && visited.insert(neighbor).second)
                        pending.push_back(neighbor);
                }
            }
        }
        return false;
    };

    std::set<Vertex> pending;
    for (Vertex v : tr.finite_vertex_handles())
        if (is_non_manifold(v)) pending.insert(v);
    std::size_t removed = 0;
    while (!pending.empty()) {
        Vertex v = *pending.begin();
        pending.erase(pending.begin());
        if (!is_non_manifold(v)) continue;
        std::vector<Cell> candidates;
        tr.finite_incident_cells(v, std::back_inserter(candidates));
        candidates.erase(std::remove_if(candidates.begin(), candidates.end(),
            [&](Cell c) { return !inside(c); }), candidates.end());
        std::stable_sort(candidates.begin(), candidates.end(), [&](Cell a, Cell b) {
            return std::abs(CGAL::to_double(tr.tetrahedron(a).volume()))
                 < std::abs(CGAL::to_double(tr.tetrahedron(b).volume()));
        });
        for (Cell c : candidates) {
            c3t3.remove_from_complex(c);
            ++removed;
            for (int i = 0; i < 4; ++i) pending.insert(c->vertex(i));
            if (!is_non_manifold(v)) break;
        }
    }

    // Replace restricted-Delaunay facets (which may include isolated sheets)
    // with exactly the boundary of the retained tetrahedra.
    for (auto f : tr.finite_facets()) {
        const bool boundary = inside(f.first) != inside(f.first->neighbor(f.second));
        if (!boundary) c3t3.remove_from_complex(f);
        else if (!c3t3.is_in_complex(f)) c3t3.add_to_complex(f, patch);
    }
    for (Vertex v : tr.finite_vertex_handles()) c3t3.set_dimension(v, -1);
    for (Cell c : c3t3.cells_in_complex())
        for (int i = 0; i < 4; ++i) {
            c3t3.set_dimension(c->vertex(i), 3);
            c3t3.set_index(c->vertex(i), subdomain);
        }
    for (auto f : c3t3.facets_in_complex())
        for (int i = 0; i < 4; ++i)
            if (i != f.second) {
                c3t3.set_dimension(f.first->vertex(i), 2);
                c3t3.set_index(f.first->vertex(i), patch);
            }
    c3t3.clear_manifold_info();
    if (c3t3.number_of_cells() == 0)
        throw std::runtime_error("Manifold cleanup removed the entire volume");
    return removed;
}

template <typename C3t3>
void write_medit(const C3t3& c3t3, const std::string& filename) {
    std::ofstream output_stream(filename);
    CGAL::IO::write_MEDIT(output_stream, c3t3,
        CGAL::parameters::all_cells(false).all_vertices(false)
    );
}

template <typename C3t3>
void write_surface(const C3t3& c3t3, const std::string& filename) {
    using Point = typename C3t3::Triangulation::Geom_traits::Point_3;
    CGAL::Surface_mesh<Point> surface;
    CGAL::facets_in_complex_3_to_triangle_mesh(c3t3, surface);
    std::ofstream output_stream(filename);
    CGAL::IO::write_OBJ(output_stream, surface, CGAL::parameters::stream_precision(17));
}

#endif // CGAL_COMMON_H
