#pragma once
// Shared data structures and constants for the PATH module.
// Converted from: Fortran PATH module COMMON blocks and parameters
// in paths.f, pathsd.f, and related files.

#include <feff/constants.hpp>
#include <feff/dimensions.hpp>
#include <feff/types.hpp>
#include <array>
#include <vector>

namespace feff::path {

// --------------------------------------------------------------------------
// Compile-time constants matching Fortran parameters
// --------------------------------------------------------------------------

// Number of energy criterion points
inline constexpr int necrit = 9;

// Number of beta angle grid points (grid from -nbeta to +nbeta)
inline constexpr int nbeta = 40;

// Heap size for path finder
inline constexpr int nx = 6000000;

// Max number of paths to save
inline constexpr int npx = 100000000;

// Number of paths to consider at one time in pathsd
inline constexpr int np1x = 600000;

// Large distance sentinel
inline constexpr float big = 1.0e3f;

// Tolerances
inline constexpr float eps3 = 1.0e-3f;   // individual leg parameters
inline constexpr float eps5 = 2.0e-5f;   // rtotal range

// --------------------------------------------------------------------------
// Shared atom data (replaces COMMON /atoms/)
// --------------------------------------------------------------------------
struct AtomData {
    // rat(3, 0:natx) -- atomic positions in single precision (path finder units)
    // Stored as rat[iat][xyz], iat from 0 to nat
    std::vector<std::array<float, 3>> rat;

    // ipot(0:natx) -- potential type index for each atom
    std::vector<int> ipot;

    // i1b(0:natx) -- first-bounce flag/degeneracy
    std::vector<int> i1b;

    void resize(int nat_plus_1) {
        rat.resize(nat_plus_1, {0.0f, 0.0f, 0.0f});
        ipot.resize(nat_plus_1, 0);
        i1b.resize(nat_plus_1, 0);
    }
};

// --------------------------------------------------------------------------
// FBeta criterion arrays (single precision, matches Fortran)
// --------------------------------------------------------------------------

// fbetac(-nbeta:nbeta, 0:nphx, necrit)
// Linearized index: [ibeta+nbeta][iph][icrit]
// Total size: (2*nbeta+1) * (nphx+1) * necrit
inline constexpr int fbetac_dim1 = 2 * nbeta + 1;
inline constexpr int fbetac_dim2 = nphx + 1;

// Index into fbetac flat array
inline int fbetac_idx(int ibeta, int iph, int icrit) {
    return (ibeta + nbeta) + fbetac_dim1 * (iph + fbetac_dim2 * icrit);
}

// fbeta(-nbeta:nbeta, 0:nphx, nex) for full energy grid
inline int fbeta_idx(int ibeta, int iph, int ie) {
    return (ibeta + nbeta) + fbetac_dim1 * (iph + fbetac_dim2 * ie);
}

} // namespace feff::path
