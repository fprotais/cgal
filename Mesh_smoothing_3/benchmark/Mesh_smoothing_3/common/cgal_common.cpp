#include "cgal_common.h"

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/boost/graph/helpers.h>
#include <CGAL/Polygon_mesh_processing/IO/polygon_mesh_io.h>
#include <CGAL/Polygon_mesh_processing/self_intersections.h>
#include <CGAL/Polygon_mesh_processing/shape_predicates.h>

bool check_input_surface(std::string input)
{
    using Kernel = CGAL::Exact_predicates_inexact_constructions_kernel;
    using Point = Kernel::Point_3;
    using Surface_mesh = CGAL::Surface_mesh<Point>;
    namespace PMP = CGAL::Polygon_mesh_processing;
    Surface_mesh surface;
    
    if (!PMP::IO::read_polygon_mesh(input, surface) || is_empty(surface) || !is_triangle_mesh(surface))
        return false;

    if (!CGAL::is_closed(surface)) return false;
    for (auto f : surface.faces())
        if (PMP::is_degenerate_triangle_face(f, surface))
            return false;
    if (PMP::does_self_intersect(surface))
        return false;
    return true;
}