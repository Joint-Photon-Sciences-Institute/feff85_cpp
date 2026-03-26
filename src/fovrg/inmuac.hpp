#pragma once
// Initialize orbital configuration and mesh for FOVRG.
// Converted from: src/FOVRG/inmuac.f

#include <feff/types.hpp>
#include <feff/dimensions.hpp>
#include "../atom/atom_types.hpp"

namespace feff::fovrg {

using feff::atom::FovrgState;

/// Initialize orbital configuration for the muffin-tin atom.
/// Calls getorb to set up orbital quantum numbers and occupations,
/// finds last tabulation points, and adds photoelectron orbital.
///
/// @param ihole  Index of core-hole orbital (0 = no hole)
/// @param xionin Ionicity
/// @param iunf   Unfolding flag for f-orbitals
/// @param ikap   Kappa quantum number for photoelectron
/// @param state  FOVRG state (modifies config, scf, orb, nuclear, lagrange, mesh)
void inmuac(int ihole, double xionin, int iunf, int ikap, FovrgState& state);

} // namespace feff::fovrg
