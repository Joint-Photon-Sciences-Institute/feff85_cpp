#pragma once
// Main driver for the complex-energy Dirac equation solver (FOVRG).
// Converted from: src/FOVRG/dfovrg.f
//
// Solves the fully relativistic Dirac equation for complex energy,
// computing regular and irregular solutions for the muffin-tin potential.

#include <feff/types.hpp>
#include <feff/dimensions.hpp>
#include "../atom/atom_types.hpp"

namespace feff::fovrg {

using feff::atom::FovrgState;

/// Main Dirac equation solver for complex energy.
///
/// @param ncycle Number of iterations with nonlocal exchange
/// @param ikap   Kappa quantum number for photoelectron
/// @param rmt    Muffin-tin radius
/// @param jlast  Last point for integration [in/out]
/// @param jri    First interstitial grid point (0-based)
/// @param p2     Complex energy [in/out, used as initial value for irregular]
/// @param dx     Loucks grid step (usually 0.05)
/// @param ri     Loucks position grid (size nrptx)
/// @param vxc    Coulomb+xc potential for total density (size nrptx)
/// @param vxcval Coulomb+xc potential for valence density (size nrptx)
/// @param dgcn   Large Dirac components for atom (nrptx x 30)
/// @param dpcn   Small Dirac components for atom (nrptx x 30)
/// @param adgc   Development coefficients for dgcn (10 x 30)
/// @param adpc   Development coefficients for dpcn (10 x 30)
/// @param xnval  Valence occupation numbers (size 30)
/// @param pu     Output: large component at muffin tin [in/out]
/// @param qu     Output: small component at muffin tin [in/out]
/// @param ps     Output: photoelectron large component (size nrptx)
/// @param qs     Output: photoelectron small component (size nrptx)
/// @param iz     Atomic number
/// @param ihole  Core-hole orbital index
/// @param xion   Ionicity
/// @param iunf   Unfolding flag for f-orbitals
/// @param irr    Solution type (<0 = regular, >0 = irregular)
/// @param ic3    Flag for c3 (spin-orbit) correction
/// @param state  FOVRG state
void dfovrg(int ncycle, int ikap, double rmt, int& jlast, int jri,
            FeffComplex& p2, double dx, const double ri[],
            FeffComplex vxc[], FeffComplex vxcval[],
            const double dgcn[][30], const double dpcn[][30],
            const double adgc[][30], const double adpc[][30],
            const double xnval[30],
            FeffComplex& pu, FeffComplex& qu,
            FeffComplex ps[], FeffComplex qs[],
            int iz, int ihole, double xion, int iunf, int irr, int ic3,
            FovrgState& state);

/// Solve Dirac equation for flat (constant) potential.
/// Given p1, q1 at r1, finds p2, q2 at r2 using exact solution
/// with spherical Bessel functions.
///
/// @param r1   First radial point
/// @param r2   Second radial point
/// @param p1   Large component at r1
/// @param q1   Small component at r1
/// @param en   Complex energy (Hartrees)
/// @param vav  Average potential (Hartrees)
/// @param ikap Kappa quantum number
/// @param p2   Output: large component at r2
/// @param q2   Output: small component at r2
void flatv(double r1, double r2, FeffComplex p1, FeffComplex q1,
           FeffComplex en, FeffComplex vav, int ikap,
           FeffComplex& p2, FeffComplex& q2);

} // namespace feff::fovrg
