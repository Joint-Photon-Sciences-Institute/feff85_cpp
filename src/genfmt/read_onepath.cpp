// Read single-path input data.
// Converted from GENFMT/read_onepath.f
//
// Stub implementation — the original Fortran used json_module.
// In the C++ version, path data is typically set programmatically
// through the OnePathInput struct.

#include "read_onepath.hpp"
#include <feff/constants.hpp>
#include <fstream>
#include <sstream>

namespace feff::genfmt {

bool read_onepath(const std::string& filename, OnePathInput& data) {
    // Initialize defaults
    data.elpty = 0.0;
    for (int i = 0; i < 3; ++i) {
        data.evec[i] = 0.0;
        data.xivec[i] = 0.0;
    }

    std::ifstream in(filename);
    if (!in.is_open()) return false;

    std::string line;

    // Read index, nleg, deg, iorder, ipol, elpty
    if (!std::getline(in, line)) return false;
    {
        std::istringstream ss(line);
        ss >> data.index >> data.nleg >> data.deg >> data.iorder
           >> data.ipol >> data.elpty;
    }

    // Read nnnn_out and json_out flags
    if (!std::getline(in, line)) return false;
    {
        std::istringstream ss(line);
        int nn, js;
        ss >> nn >> js;
        data.nnnn_out = (nn > 0);
        data.json_out = (js > 0);
    }

    // Read atom coordinates (in Angstrom) and ipot
    for (int iat = 1; iat <= data.nleg; ++iat) {
        if (!std::getline(in, line)) return false;
        std::istringstream ss(line);
        double x, y, z;
        int ip;
        ss >> x >> y >> z >> ip;
        // Convert to code units (Bohr)
        data.rat[0][iat] = x / bohr;
        data.rat[1][iat] = y / bohr;
        data.rat[2][iat] = z / bohr;
        data.ipot[iat] = ip;
    }

    // Read evec and xivec
    if (std::getline(in, line)) {
        std::istringstream ss(line);
        ss >> data.evec[0] >> data.evec[1] >> data.evec[2];
    }
    if (std::getline(in, line)) {
        std::istringstream ss(line);
        ss >> data.xivec[0] >> data.xivec[1] >> data.xivec[2];
    }

    return true;
}

bool read_geometry(const std::string& filename, OnePathInput& data) {
    std::ifstream in(filename);
    if (!in.is_open()) return false;

    std::string line;

    // Read nleg
    if (!std::getline(in, line)) return false;
    {
        std::istringstream ss(line);
        ss >> data.nleg;
    }

    // Read atom coordinates
    for (int iat = 1; iat <= data.nleg; ++iat) {
        if (!std::getline(in, line)) return false;
        std::istringstream ss(line);
        double x, y, z;
        int ip;
        ss >> x >> y >> z >> ip;
        data.rat[0][iat] = x / bohr;
        data.rat[1][iat] = y / bohr;
        data.rat[2][iat] = z / bohr;
        data.ipot[iat] = ip;
    }

    return true;
}

} // namespace feff::genfmt
