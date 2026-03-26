// Heap operations for the path finder.
// Converted from: src/PATH/heap.f
//
// IMPORTANT: Arrays use 1-based indexing to match the Fortran exactly.
// h[1..n] is the heap, ih[1..n] are associated indices.
// Element 0 is unused (caller allocates n+1 elements).

#include "heap.hpp"
#include <utility>

namespace feff::path {

void hup(float h[], int ih[], int n) {
    // New element is at position n; bubble it up.
    int i = n;

    for (;;) {
        int j = i / 2;
        // If no parent, we're at the top -- done
        if (j == 0) return;
        if (h[i] < h[j]) {
            std::swap(h[i], h[j]);
            std::swap(ih[i], ih[j]);
            i = j;
        } else {
            return;
        }
    }
}

void hdown(float h[], int ih[], int n) {
    // Element at position 1 has been replaced; bubble it down.
    int i = 1;

    for (;;) {
        int j = 2 * i;
        int k = j + 1;

        // If j > n, element is at bottom -- done
        if (j > n) return;
        // Handle case where element has only one child
        if (k > n) k = j;

        // j = index of smallest child
        if (h[j] > h[k]) j = k;

        if (h[i] > h[j]) {
            std::swap(h[i], h[j]);
            std::swap(ih[i], ih[j]);
            i = j;
        } else {
            return;
        }
    }
}

} // namespace feff::path
