#pragma once
// Map edge name to core-hole index.
// Converted from: src/RDINP/setedg.f

#include <string>

namespace feff::rdinp {

/// Convert edge name (e.g. "K", "L1", "L3") or numeric string to ihole index.
/// Throws std::runtime_error if the edge name is not recognized.
/// Replaces Fortran subroutine setedg(a2, ihole).
void setedg(const std::string& edge_name, int& ihole);

} // namespace feff::rdinp
