#pragma once

// GENFMT initialization (reads phase.pad, sets up arrays).
// Converted from GENFMT/genfmt_prep.f

#include "genfmt_data.hpp"
#include <string>

namespace feff::genfmt {

/// Prepare all data needed for the genfmt calculation.
/// Reads phase.pad, calls rdxsph (external), setkap, snlm, and sets up
/// arrays for xk, ck, ckmag, xkr.
///
/// phpad: path to phase.pad file
/// ispin: spin flag
///
/// Output: fills PhaseData, SpinPhaseData, NlmData, and various arrays.
/// Also sets: nsp, ll, npath, ntotal, nused, xportx
struct GenfmtPrepResult {
    // Energy grid info
    int ne = 0;
    int ne1 = 0;
    int ne3 = 0;
    int ik0 = 0;
    int ixc = 0;

    // Potential info
    int npot = 0;
    int ihole = 0;
    double rnrmav = 0.0;
    double xmu = 0.0;
    double edge = 0.0;
    double rs = 0.0;
    double vint = 0.0;

    // Angular momentum
    int kinit = 0;
    int linit = 0;
    int ilinit = 0;
    int lmaxp1 = 0;

    // Spin
    int nsp = 1;
    int ll = 0;

    // Path counters
    int npath = 0;
    int ntotal = 0;
    int nused = 0;
    double xportx = -1.0;

    // Arrays
    double xk[nex]{};
    double ckmag[nex]{};
    double xkr[nex]{};
    FeffComplex ck[nex]{};
};

/// Initialize genfmt data by reading phase.pad and setting up arrays.
///
/// phpad: path to phase.pad
/// ispin: spin flag
/// pd: phase data (output)
/// spd: spin-resolved phase data (output)
/// nlm: Legendre normalization factors (output)
/// prep: preparation result with all computed arrays (output)
void genfmt_prep(const std::string& phpad, int ispin,
                 PhaseData& pd, SpinPhaseData& spd, NlmData& nlm,
                 GenfmtPrepResult& prep);

} // namespace feff::genfmt
