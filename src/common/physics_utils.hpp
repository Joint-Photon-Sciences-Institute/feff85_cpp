#pragma once
// Physics utility functions used across multiple FEFF modules.
// Converted from: src/COMMON/getxk.f, pijump.f, setkap.f, setgam.f, iniptz.f

#include <feff/types.hpp>

namespace feff::common {

/// Calculate wavenumber k from energy E (in Hartrees).
///   k = sqrt(2*E) for E > 0 (above edge)
///   k = -sqrt(-2*E) for E < 0 (below edge)
/// Replaces Fortran getxk(e).
double getxk(double e);

/// Remove 2π jumps in phase to maintain continuity.
/// ph may be adjusted by multiples of 2π. old is the previous phase value.
/// Replaces Fortran pijump(ph, old).
void pijump(double& ph, double old);

/// Set initial state angular momentum l and kappa quantum number from core-hole index.
///   ihole: 1=K(1s), 2=LI(2s), 3=LII(2p1/2), 4=LIII(2p3/2), ...
/// Replaces Fortran setkap(ihole, kinit, linit).
void setkap(int ihole, int& kinit, int& linit);

/// Set core-hole lifetime broadening (Gamma_ch) in eV from Rahkonen/Krause tables.
/// Replaces Fortran setgam(iz, ihole, gamach).
void setgam(int iz, int ihole, double& gamach);

/// Initialize polarization tensor for X-ray spectroscopy.
/// ptz is a 3x3 complex matrix indexed as ptz[i+1][j+1] for i,j in {-1,0,1}.
///   iptz: polarization type (1-10, where 10=orientation average)
///   modus: 1=spherical coordinates, 2=Cartesian coordinates
/// Replaces Fortran iniptz(ptz, iptz, modus).
void iniptz(FeffComplex ptz[3][3], int iptz, int modus);

} // namespace feff::common
