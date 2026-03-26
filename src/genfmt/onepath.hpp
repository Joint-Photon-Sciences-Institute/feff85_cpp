#pragma once

// Single-path F-matrix calculation.
// Converted from GENFMT/onepath.f
//
// Computes the F-matrix for a single scattering path, returning
// the columns of the feffNNNN.dat file.

#include "genfmt_data.hpp"
#include <string>

namespace feff::genfmt {

/// Result columns from a single-path calculation (feffNNNN.dat columns).
struct OnePathResult {
    int ne1 = 0;            // Number of k-grid points
    double col1[nex]{};     // k-grid
    double col2[nex]{};     // central atom phase shifts
    double col3[nex]{};     // magnitude of F_eff
    double col4[nex]{};     // phase of F_eff
    double col5[nex]{};     // reduction factor
    double col6[nex]{};     // mean free path
    double col7[nex]{};     // real part of complex momentum

    // Geometry outputs
    double ri[legtot]{};
    double beta[legtot + 1]{};
    double eta[legtot + 2]{};
};

/// Potential-related output from onepath.
struct OnePathPotInfo {
    std::string cxc;        // Exchange-correlation model description
    double rs = 0.0;        // Interstitial radius
    double vint = 0.0;      // Interstitial potential
    double xmu = 0.0;       // Fermi energy
    double edge = 0.0;      // Threshold relative to atomic value
    double xkf = 0.0;       // k at Fermi energy
    double rnrmav = 0.0;    // Average Norman radius
    double gamach = 0.0;    // Core-hole lifetime broadening
    std::string versn;      // Version string
};

/// Compute F-matrix for a single path.
///
/// phpad: path to phase.pad
/// index: path index
/// nleg: number of legs
/// deg: path degeneracy
/// iorder: approximation order
/// ipot: potential indices for each atom
/// rat: atom positions [3][legtot+2] (in Angstrom, converted internally to Bohr)
/// iz: atomic numbers
/// ipol: polarization flag
/// evec: polarization vector
/// elpty: ellipticity
/// xivec: direction of travel
/// write_feffdat: write feffNNNN.dat file
/// write_xdi: write feffNNNN.xdi file
/// verbose: write screen messages
///
/// Output: result contains the 7 columns and geometry data.
///         potinfo contains potential-related information.
void onepath(const std::string& phpad, int index, int nleg, double deg,
             int iorder, int ipot[], double rat[][legtot + 2],
             int iz[],
             int ipol, const double evec[3], double elpty,
             const double xivec[3],
             bool write_feffdat, bool write_xdi, bool verbose,
             OnePathResult& result, OnePathPotInfo& potinfo);

} // namespace feff::genfmt
