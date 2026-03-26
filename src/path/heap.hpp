#pragma once
// Heap operations for the path finder.
// Converted from: src/PATH/heap.f

namespace feff::path {

/// Bubble last element up to its proper location in the min-heap.
/// h[0..n-1] is the heap, ih[0..n-1] are associated indices.
/// n is the current heap size (1-based count; arrays are 0-based internally
/// but we keep 1-based indexing to match Fortran exactly).
void hup(float h[], int ih[], int n);

/// Bubble element at position 1 down to its proper location.
/// Element at h[1] has been replaced; restore heap property.
void hdown(float h[], int ih[], int n);

} // namespace feff::path
