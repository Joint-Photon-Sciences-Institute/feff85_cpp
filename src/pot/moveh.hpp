#pragma once
// Adjust hydrogen atom positions for better MT geometry.
// Converted from src/POT/moveh.f
//
// Increases bond lengths for hydrogen atoms to avoid screwing up
// muffin-tin geometry.  The maximum distance is set empirically
// from H2O and GeH4 calculations.

#include <feff/dimensions.hpp>

namespace feff::pot {

/// Move hydrogen atoms outward to improve MT sphere geometry.
///
/// @param nat    Number of atoms in the cluster
/// @param iphat  Potential type for each atom [natx] (0-based atom index)
/// @param iz     Atomic number per potential type [nphx+1] (0-based pot index)
/// @param rath   Atomic coordinates [3][natx] (column-major, modified in place)
void moveh(int nat, const int* iphat, const int* iz, double* rath);

} // namespace feff::pot
