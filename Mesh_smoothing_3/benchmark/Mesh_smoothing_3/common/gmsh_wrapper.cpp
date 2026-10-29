#include "gmsh_wrapper.h"

#include <gmsh.h>

#include <algorithm>
#include <cmath>
#include <filesystem>
#include <limits>
#include <mutex>
#include <stdexcept>
#include <vector>

namespace {

// Gmsh uses global state. Serialize calls to this helper, and preserve any
// session/models that were created by the caller.
std::mutex gmsh_mutex;

class Gmsh_session {
public:
    Gmsh_session() : owns_api(!gmsh::isInitialized()) {
        if (owns_api) gmsh::initialize(0, nullptr, false);
    }

    void read_mesh(const std::string& filename) {
        gmsh::option::getNumber("General.Terminal", previous_terminal);
        gmsh::option::getNumber("General.AbortOnError", previous_abort_on_error);
        options_saved = true;
        gmsh::option::setNumber("General.Terminal", 0);
        gmsh::option::setNumber("General.AbortOnError", 2);
        gmsh::model::getCurrent(previous_model);

        std::vector<std::string> models;
        gmsh::model::list(models);
        model_name = "benchmark_gmsh_metrics";
        while (std::find(models.begin(), models.end(), model_name) != models.end())
            model_name += '_';
        gmsh::model::add(model_name);
        model_created = true;
        gmsh::merge(filename);
    }

    ~Gmsh_session() {
        try {
            if (owns_api) {
                gmsh::finalize();
            } else {
                if (model_created) {
                    gmsh::model::setCurrent(model_name);
                    gmsh::model::remove();
                }
                if (!previous_model.empty()) gmsh::model::setCurrent(previous_model);
                if (options_saved) {
                    gmsh::option::setNumber("General.Terminal", previous_terminal);
                    gmsh::option::setNumber("General.AbortOnError", previous_abort_on_error);
                }
            }
        } catch (...) {
            // Cleanup must not replace an exception from reading/evaluating the mesh.
        }
    }

    Gmsh_session(const Gmsh_session&) = delete;
    Gmsh_session& operator=(const Gmsh_session&) = delete;

private:
    bool owns_api;
    bool options_saved = false;
    bool model_created = false;
    double previous_terminal = 0;
    double previous_abort_on_error = 0;
    std::string previous_model;
    std::string model_name;
};

} // namespace

Gmsh_metrics compute_gmsh_metrics(const std::string& filename) {
    const std::lock_guard<std::mutex> lock(gmsh_mutex);
    const std::string error_prefix = "Cannot compute Gmsh SICN metrics for '" + filename + "': ";
    try {
        Gmsh_session session;
        session.read_mesh(filename);

        std::vector<int> element_types;
        std::vector<std::vector<std::size_t>> element_tags, node_tags;
        gmsh::model::mesh::getElements(element_types, element_tags, node_tags, 3);

        Gmsh_metrics metrics{std::numeric_limits<double>::max(),
                             std::numeric_limits<double>::lowest(), 0};
        long double sum = 0;
        std::size_t count = 0;
        for (const auto& tags : element_tags) {
            if (tags.empty()) continue;
            std::vector<double> qualities;
            gmsh::model::mesh::getElementQualities(tags, qualities, "minSICN");
            if (qualities.size() != tags.size())
                throw std::runtime_error("Gmsh returned an incomplete set of element qualities");
            for (double quality : qualities) {
                if (!std::isfinite(quality))
                    throw std::runtime_error("Gmsh returned a non-finite element quality");
                metrics.SICN_min = (std::min)(metrics.SICN_min, quality);
                metrics.SICN_max = (std::max)(metrics.SICN_max, quality);
                sum += quality;
                ++count;
            }
        }
        if (count == 0) throw std::runtime_error("mesh contains no 3D elements");
        metrics.SICN_avg = static_cast<double>(sum / count);
        return metrics;
    } catch (const std::exception& error) {
        throw std::runtime_error(error_prefix + error.what());
    } catch (const std::string& error) {
        // Some versions of the native Gmsh API throw strings.
        throw std::runtime_error(error_prefix + error);
    } catch (int error) {
        throw std::runtime_error(error_prefix + "Gmsh error code " + std::to_string(error));
    } catch (...) {
        throw std::runtime_error(error_prefix + "unknown Gmsh error");
    }
}

void vtk2medit(const std::string& vtk_filename, const std::string& medit_filename) {
    const std::lock_guard<std::mutex> lock(gmsh_mutex);
    const std::string error_prefix = "Cannot convert VTK to Medit for '" + vtk_filename + "': ";
    try {
        // Gmsh selects its writer from the output filename's extension.
        if (std::filesystem::path(medit_filename).extension() != ".mesh")
            throw std::runtime_error("Medit output filename must end in .mesh");
        Gmsh_session session;
        session.read_mesh(vtk_filename);

        std::vector<int> element_types;
        gmsh::model::mesh::getElementTypes(element_types);
        if (element_types.empty()) throw std::runtime_error("input mesh contains no elements");
        for (int type : element_types) {
            // Gmsh's Medit writer supports only these linear element types.
            // In particular, it silently skips prisms, pyramids and higher-order cells.
            if (type != 1 && type != 2 && type != 3 && type != 4 && type != 5)
                throw std::runtime_error("Gmsh's Medit writer does not support element type "
                                         + std::to_string(type));
        }
        gmsh::write(medit_filename);
    } catch (const std::exception& error) {
        throw std::runtime_error(error_prefix + error.what());
    } catch (const std::string& error) {
        throw std::runtime_error(error_prefix + error);
    } catch (int error) {
        throw std::runtime_error(error_prefix + "Gmsh error code " + std::to_string(error));
    } catch (...) {
        throw std::runtime_error(error_prefix + "unknown Gmsh error");
    }
}
