#pragma once
// Heap-sort routines that return sorted index arrays.
// Converted from: src/PATH/sortix.f
// Following algorithm in Knuth, Vol 3, pp 146-7.

namespace feff::path {

/// Sort by rearranging indices; keys are float (real) numbers.
/// r[1..n] are the keys (1-based). index[1..n] is returned such that
/// r[index[1]] <= r[index[2]] <= ... <= r[index[n]].
/// r is NOT modified.
void sortir(int n, int index[], const float r[]);

/// Sort by rearranging indices; keys are integers.
/// k[1..n] are the keys (1-based). index[1..n] returned sorted.
void sortii(int n, int index[], const int k[]);

/// Sort by rearranging indices; keys are double precision numbers.
/// r[1..n] are the keys (1-based). index[1..n] returned sorted.
void sortid(int n, int index[], const double r[]);

} // namespace feff::path
