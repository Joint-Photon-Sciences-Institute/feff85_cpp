// Heap-sort routines returning sorted index arrays.
// Converted from: src/PATH/sortix.f
// Knuth, The Art of Computer Programming, Vol 3, pp 146-7.
//
// IMPORTANT: All arrays are 1-based to match Fortran exactly.
// Caller must allocate index[n+1] and r[n+1] with element 0 unused.

#include "sortix.hpp"

namespace feff::path {

void sortir(int n, int index[], const float r[]) {
    // Initialize index array (1-based)
    for (int i = 1; i <= n; ++i) {
        index[i] = i;
    }
    if (n == 1) return;

    // H1: initialize
    int l = n / 2 + 1;
    int ir = n;

    for (;;) {
        int irr;
        float rr;

        // H2: Decrease l or ir
        if (l > 1) {
            --l;
            irr = index[l];
            rr = r[irr];
        } else {
            irr = index[ir];
            rr = r[irr];
            index[ir] = index[1];
            --ir;
            if (ir == 1) {
                index[1] = irr;
                return;
            }
        }

        // H3: Prepare for sift-up
        int j = l;

        for (;;) {
            // H4: Advance downward
            int i = j;
            j = 2 * j;
            if (j == ir) {
                // H6: Son larger than rr?
                if (rr >= r[index[j]]) {
                    // H8
                    index[i] = irr;
                    break;
                }
                // H7: Move son up
                index[i] = index[j];
                // continue sifting (but j will be > ir next iteration)
                continue;
            }
            if (j > ir) {
                // H8: Store rr
                index[i] = irr;
                break;
            }

            // H5: Find larger son of i
            if (r[index[j]] < r[index[j + 1]]) ++j;

            // H6: Son larger than rr?
            if (rr >= r[index[j]]) {
                // H8
                index[i] = irr;
                break;
            }

            // H7: Move son up
            index[i] = index[j];
        }
    }
}

void sortii(int n, int index[], const int k[]) {
    for (int i = 1; i <= n; ++i) {
        index[i] = i;
    }
    if (n == 1) return;

    int l = n / 2 + 1;
    int ir = n;

    for (;;) {
        int irr, kk;

        if (l > 1) {
            --l;
            irr = index[l];
            kk = k[irr];
        } else {
            irr = index[ir];
            kk = k[irr];
            index[ir] = index[1];
            --ir;
            if (ir == 1) {
                index[1] = irr;
                return;
            }
        }

        int j = l;

        for (;;) {
            int i = j;
            j = 2 * j;
            if (j == ir) {
                if (kk >= k[index[j]]) {
                    index[i] = irr;
                    break;
                }
                index[i] = index[j];
                continue;
            }
            if (j > ir) {
                index[i] = irr;
                break;
            }
            if (k[index[j]] < k[index[j + 1]]) ++j;
            if (kk >= k[index[j]]) {
                index[i] = irr;
                break;
            }
            index[i] = index[j];
        }
    }
}

void sortid(int n, int index[], const double r[]) {
    for (int i = 1; i <= n; ++i) {
        index[i] = i;
    }
    if (n == 1) return;

    int l = n / 2 + 1;
    int ir = n;

    for (;;) {
        int irr;
        double rr;

        if (l > 1) {
            --l;
            irr = index[l];
            rr = r[irr];
        } else {
            irr = index[ir];
            rr = r[irr];
            index[ir] = index[1];
            --ir;
            if (ir == 1) {
                index[1] = irr;
                return;
            }
        }

        int j = l;

        for (;;) {
            int i = j;
            j = 2 * j;
            if (j == ir) {
                if (rr >= r[index[j]]) {
                    index[i] = irr;
                    break;
                }
                index[i] = index[j];
                continue;
            }
            if (j > ir) {
                index[i] = irr;
                break;
            }
            if (r[index[j]] < r[index[j + 1]]) ++j;
            if (rr >= r[index[j]]) {
                index[i] = irr;
                break;
            }
            index[i] = index[j];
        }
    }
}

} // namespace feff::path
