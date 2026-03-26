#pragma once

// Read genfmt configuration (replacement for JSON/input file reading).
// Converted from GENFMT/regenf.f
//
// In the C++ version, this reads configuration from a simple settings
// structure rather than JSON files.

#include <feff/types.hpp>

namespace feff::genfmt {

/// Configuration for the genfmt module.
struct GenfmtConfig {
    int mfeff = 1;       // Flag for EXAFS calculation
    int ipr5 = 0;        // Print level for feff.dat output
    int iorder = 2;      // Order of approximation in f-matrix expansion
    bool wnstar = false; // Write nstar.dat flag
    double critcw = 4.0; // Curved wave chi amplitude ratio threshold (%)

    // Polarization parameters (from global.dat equivalent)
    int ipol = 0;        // Polarization type (0=average)
    int ispin = 0;       // Spin flag
    int le2 = 0;         // Multipole moment selector
    double angks = 0.0;  // Angle between k-vector and spin-vector
    double elpty = 0.0;  // Ellipticity
    double evec[3] = {}; // Polarization vector
    double xivec[3] = {};// Direction of travel
    FeffComplex ptz[3][3] = {}; // Polarization tensor (indexed [i+1][j+1])
};

/// Read genfmt configuration.
/// In the original code, this read from genfmt.json and global.json.
/// Here it populates the GenfmtConfig structure.
void regenf(GenfmtConfig& config);

} // namespace feff::genfmt
