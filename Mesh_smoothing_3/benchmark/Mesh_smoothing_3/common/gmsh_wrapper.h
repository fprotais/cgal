#ifndef GMSH_METRICS_H
#define GMSH_METRICS_H

#include <iostream>
#include <string>
#include <sstream>


struct Gmsh_metrics {
    double SICN_min;
    double SICN_max;
    double SICN_avg;
};

Gmsh_metrics compute_gmsh_metrics(const std::string& filename);

void vtk2medit(const std::string& vtk_filename, const std::string& medit_filename);

#endif // GMSH_METRICS_H

