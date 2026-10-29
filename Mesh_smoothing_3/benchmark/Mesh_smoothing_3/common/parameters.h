#ifndef PARAMETERS_H
#define PARAMETERS_H

#include <iostream>
#include <string>
#include <sstream>
#include <type_traits>




struct Input_parameters {
    double boundary_weight = 1.;
    bool exact_checks = false;
    unsigned grid_size = 100;
    double bbox_target_sizing = 0.05;
    unsigned octree_cube_cells_per_dim = 10;
    unsigned octree_nb_levels = 2;
    double offset = 0.01;
    bool verbose = false;

    friend std::ostream& operator<<(std::ostream& os, const Input_parameters& params) {
        os << "    Boundary Weight: " << params.boundary_weight << "\n";
        os << "    Exact Checks: " << (params.exact_checks ? "true" : "false") << "\n";
        os << "    Grid Size: " << params.grid_size << "\n";
        os << "    BBox Target Sizing: " << params.bbox_target_sizing << "\n";
        os << "    Octree Cube Cells Per Dim: " << params.octree_cube_cells_per_dim << "\n";
        os << "    Octree Nb Levels: " << params.octree_nb_levels << "\n";
        os << "    Offset: " << params.offset << "\n";
        os << "    Verbose: " << (params.verbose ? "true" : "false") << "\n";
        return os;
    }

    // Read name value or name = value lines; blank lines and # comments are allowed.
    // Invalid entries are reported and skipped, keeping the field's current value.
    friend std::istream& operator>>(std::istream& is, Input_parameters& params) {
        std::string line;
        std::size_t line_number = 0;
        while (std::getline(is, line)) {
            ++line_number;
            const auto comment = line.find('#');
            if (comment != std::string::npos) line.erase(comment);
            if (line.find_first_not_of(" \t\r") == std::string::npos) continue;

            const auto equals = line.find('=');
            if (equals != std::string::npos) line[equals] = ' ';
            const auto warn = [&] {
                std::cerr << "Warning: ignoring invalid parameter at line "
                          << line_number << ": " << line << '\n';
            };
            std::istringstream entry(line);
            std::string key, value, extra;
            if (!(entry >> key >> value) || (entry >> extra)) {
                warn();
                continue;
            }

            // Parse and validate before assigning, using the field's own type.
            const auto read_value = [&](auto& field, bool strictly_positive = false) {
                using Value = std::decay_t<decltype(field)>;
                Value parsed{};
                if constexpr (std::is_same_v<Value, bool>) {
                    if (value == "true" || value == "1") parsed = true;
                    else if (value == "false" || value == "0") parsed = false;
                    else return false;
                } else {
                    if constexpr (std::is_unsigned_v<Value>)
                        if (value.front() == '-') return false;
                    std::istringstream input(value);
                    if (!(input >> parsed)) return false;
                    input >> std::ws;
                    if (!input.eof()) return false;
                    if (!std::isfinite(parsed) || parsed < 0 || (strictly_positive && parsed == 0))
                        return false;
                }
                field = parsed;
                return true;
            };

            // Add one branch here for each supported parameter.
            bool valid = false;
            if (key == "boundary_weight") valid = read_value(params.boundary_weight);
            else if (key == "exact_checks") valid = read_value(params.exact_checks);
            else if (key == "grid_size") valid = read_value(params.grid_size, true);
            else if (key == "bbox_target_sizing") valid = read_value(params.bbox_target_sizing, true);
            else if (key == "octree_cube_cells_per_dim") valid = read_value(params.octree_cube_cells_per_dim, true);
            else if (key == "octree_nb_levels") valid = read_value(params.octree_nb_levels);
            else if (key == "offset") valid = read_value(params.offset, true);
            else if (key == "verbose") valid = read_value(params.verbose);
            if (!valid) {
                warn();
            }
        }
        if (is.eof() && !is.bad())
            is.clear(is.rdstate() & ~std::ios::failbit);
        return is;
    }

    static std::string csv_header() {
        return "Boundary Weight,Exact Checks,Grid Size,BBox Target Sizing,Octree Cube Cells Per Dim,Octree Nb Levels,Offset,Verbose";
    }

    std::string csv_metric() const {
        std::ostringstream oss;
        oss << boundary_weight << ",";
        oss << (exact_checks ? "true" : "false") << ",";
        oss << grid_size << ",";
        oss << bbox_target_sizing << ",";
        oss << octree_cube_cells_per_dim << ",";
        oss << octree_nb_levels << ",";
        oss << offset << ",";
        oss << (verbose ? "true" : "false");
        return oss.str();
    }

    static std::string default_file() {
        std::ostringstream oss;
        oss << "boundary_weight " << 1.0 << "\n";
        oss << "exact_checks " << true << "\n";
        oss << "grid_size " << 10 << "\n";
        oss << "bbox_target_sizing " << 0.1 << "\n";
        oss << "octree_cube_cells_per_dim " << 10 << "\n";
        oss << "octree_nb_levels " << 2 << "\n";
        oss << "offset " << 0.01 << "\n";
        oss << "verbose " << false << "\n";
        return oss.str();
    }
};


#endif // PARAMETERS_H
