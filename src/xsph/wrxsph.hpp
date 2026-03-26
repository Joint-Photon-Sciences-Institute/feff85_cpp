#pragma once
// Write phase.pad output in PAD format.
// Converted from src/XSPH/wrxsph.f

#include <feff/types.hpp>
#include <feff/dimensions.hpp>
#include <string>

namespace feff::xsph {

/// Write phase.pad file for paths and genfmt modules.
///
/// @param phpad    Output filename (e.g., "phase.pad")
/// @param nsp     Number of spin channels (1 or 2)
/// @param ne      Total energy points
/// @param ne1     Horizontal grid points
/// @param ne3     Auxiliary horizontal points
/// @param nph     Number of unique potentials
/// @param ihole   Core-hole orbital index
/// @param rnrmav  Average Norman radius
/// @param xmu     Fermi energy (Hartrees)
/// @param edge    X-ray frequency at Fermi level (Hartrees)
/// @param ik0     Grid index at Fermi level
/// @param ixc     Potential model
/// @param rs      Density parameter rs
/// @param vint    Muffin-tin zero (Hartrees)
/// @param em      Complex energy grid [nex]
/// @param eref    Energy reference [nex * nspx]
/// @param lmax    Max angular momentum per potential [nphx+1]
/// @param iz      Atomic numbers [nphx+1]
/// @param potlbl  Potential labels [nphx+1]
/// @param ph      Phase shifts [nex * (2*ltot+1) * nspx * (nphx+1)]
/// @param rkk     Multipole matrix elements [nex * 8 * nspx]
void wrxsph(const std::string& phpad,
            int nsp, int ne, int ne1, int ne3, int nph, int ihole,
            double rnrmav, double xmu, double edge,
            int ik0, int ixc, double rs, double vint,
            const FeffComplex em[], const FeffComplex* eref,
            const int lmax[], const int iz[], const char potlbl[][7],
            const FeffComplex* ph, const FeffComplex* rkk);

} // namespace feff::xsph
