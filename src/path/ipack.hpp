#pragma once
// Integer packing/unpacking for path atom indices.
// Converted from: src/PATH/ipack.f
//
// Packs n (<=8) atom indices (each 0..1289) into 3 integers,
// using base-1290 encoding.

#include <array>

namespace feff::path {

/// Pack n atom indices into 3 integers.
/// ipat[0..n-1] are the atom indices (each 0..1289).
/// iout[0..2] receives the packed result.
void ipack(int iout[3], int n, const int ipat[]);

/// Unpack 3 integers back to n atom indices.
/// On input, n is the max number to retrieve (must be <= 8).
/// On output, n is the actual count unpacked.
/// ipat[0..n-1] receives the atom indices.
void upack(const int iout[3], int& n, int ipat[]);

} // namespace feff::path
