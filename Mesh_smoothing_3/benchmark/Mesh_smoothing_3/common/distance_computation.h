#ifndef DISTANCE_COMPUTATION_H
#define DISTANCE_COMPUTATION_H

#include <string>

struct Distance_metrics {
    double hausdorff_distance;
    double l2_distance;
    double normal_deviation;
    double features_hausdorff_distance;
    double features_l2_distance;
    double corner_hausdorff_distance;
};


// mainly written by codex, there might be mistake still, need checking
Distance_metrics compute_distance_metrics(const std::string& reference_surface, const std::string& generated_surface);

Distance_metrics compute_offset_distance_metrics(const std::string& reference_surface,
                                                const std::string& generated_surface,
                                                double offset);

#endif // DISTANCE_COMPUTATION_H
