#ifndef STATISTICS_H
#define STATISTICS_H

#include <iostream>
#include <string>
#include <sstream>
#include <cmath>
#include <chrono>
#include "distance_computation.h"
#include "gmsh_wrapper.h"


struct Timer {
    std::chrono::steady_clock::time_point start_time;
    Timer() {
        start_time = std::chrono::steady_clock::now();
    }

    double pop() {
        auto end_time = std::chrono::steady_clock::now();
        return std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time).count() / 1000.0;
    }
};



struct Mesh_metrics {
    double computation_time;

    unsigned nb_inverted_elements;
    Gmsh_metrics gmsh_metrics;
    Distance_metrics distance_metrics;

    friend std::ostream& operator<<(std::ostream& os, const Mesh_metrics& metric) {
        os << "    Computation Time: " << metric.computation_time << "\n";
        os << "    SICN Min: " << metric.gmsh_metrics.SICN_min << "\n";
        os << "    SICN Max: " << metric.gmsh_metrics.SICN_max << "\n";
        os << "    SICN Avg: " << metric.gmsh_metrics.SICN_avg << "\n";
        os << "    Hausdorff Distance: " << metric.distance_metrics.hausdorff_distance << "\n";
        os << "    L2 Distance: " << metric.distance_metrics.l2_distance << "\n";
        os << "    Normal Deviation: " << metric.distance_metrics.normal_deviation << "\n";
        os << "    Features Hausdorff Distance: " << metric.distance_metrics.features_hausdorff_distance << "\n";
        os << "    Features L2 Distance: " << metric.distance_metrics.features_l2_distance << "\n";
        os << "    Corner Hausdorff Distance: " << metric.distance_metrics.corner_hausdorff_distance << "\n";
        os << "    Inverted Elements: " << metric.nb_inverted_elements << "\n";
        return os;
    }

    static std::string csv_header(std::string prefix) {
        return prefix + "Computation Time,"
             + prefix + "SICN Min,"
             + prefix + "SICN Max,"
             + prefix + "SICN Avg,"
             + prefix + "Hausdorff Distance,"
             + prefix + "L2 Distance,"
             + prefix + "Normal Deviation,"
             + prefix + "Features Hausdorff Distance,"
             + prefix + "Features L2 Distance,"
             + prefix + "Corner Hausdorff Distance,"
             + prefix + "Inverted Elements";
    }

    std::string csv_metric() const {
        std::ostringstream oss;
        oss << computation_time << ",";
        oss << gmsh_metrics.SICN_min << ",";
        oss << gmsh_metrics.SICN_max << ",";
        oss << gmsh_metrics.SICN_avg << ",";
        oss << distance_metrics.hausdorff_distance << ",";
        oss << distance_metrics.l2_distance << ",";
        oss << distance_metrics.normal_deviation << ",";
        oss << distance_metrics.features_hausdorff_distance << ",";
        oss << distance_metrics.features_l2_distance << ",";
        oss << distance_metrics.corner_hausdorff_distance << ",";
        oss << nb_inverted_elements;
        return oss.str();
    }
};

struct Result_and_statistics {
    std::string run_type;
    std::string case_name;
    unsigned nb_vertices;
    unsigned nb_faces;
    unsigned nb_tetrahedra;
    double time_preprocessing;
    Mesh_metrics init_metrics;
    Mesh_metrics smoothed_metrics;

    friend std::ostream& operator<<(std::ostream& os, const Result_and_statistics& result) {
        os << "Run Type: " << result.run_type << "\n";
        os << "Case Name: " << result.case_name << "\n";
        os << "Vertices: " << result.nb_vertices << "\n";
        os << "Faces: " << result.nb_faces << "\n";
        os << "Tetrahedra: " << result.nb_tetrahedra << "\n";
        os << "Preprocessing Time: " << result.time_preprocessing << "\n";
        os << "Initial Metrics: \n" << result.init_metrics << "\n";
        os << "Smoothed Metrics: \n" << result.smoothed_metrics << "\n";
        return os;
    }

    static std::string csv_header() {
        std::ostringstream oss;
        oss << "Run Type,Case Name,Vertices,Faces,Tetrahedra,Preprocessing Time,";
        oss << Mesh_metrics::csv_header("Initial ") << ",";
        oss << Mesh_metrics::csv_header("Smoothed ") << ",";
        return oss.str();
    }

    std::string csv_metric() const {
        std::ostringstream oss;
        oss << run_type << ",";
        oss << case_name << ",";
        oss << nb_vertices << ",";
        oss << nb_faces << ",";
        oss << nb_tetrahedra << ",";
        oss << time_preprocessing << ",";
        oss << init_metrics.csv_metric() << ",";
        oss << smoothed_metrics.csv_metric() << ",";
        return oss.str();
    }
};


#endif // STATISTICS_H
