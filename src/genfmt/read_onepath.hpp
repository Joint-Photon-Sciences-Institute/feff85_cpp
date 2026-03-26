#pragma once

// Read single-path input data.
// Converted from GENFMT/read_onepath.f
//
// The original Fortran used json_module to read feffpath.json.
// In C++, we provide both a struct-based interface and optional
// JSON reading capability.

#include <feff/dimensions.hpp>
#include <string>

namespace feff::genfmt {

/// Data for a single path, equivalent to what json_read_onepath reads.
struct OnePathInput {
    int index = 0;
    int nleg = 0;
    double deg = 0.0;
    int iorder = 2;
    int ipol = 0;
    double elpty = 0.0;
    double evec[3] = {};
    double xivec[3] = {};
    bool nnnn_out = false;
    bool json_out = false;

    // Atom coordinates in Angstrom
    double rat[3][legtot + 2]{};
    // Potential indices
    int ipot[legtot + 1]{};
};

/// Read single-path geometry from a file.
/// In the Fortran code this read from feffpath.json.
/// The C++ version reads from a simple text format or can be
/// populated programmatically.
///
/// filename: input file path
/// data: output path data
/// Returns true on success.
bool read_onepath(const std::string& filename, OnePathInput& data);

/// Read only the geometry portion (nleg, rat, ipot).
bool read_geometry(const std::string& filename, OnePathInput& data);

} // namespace feff::genfmt
