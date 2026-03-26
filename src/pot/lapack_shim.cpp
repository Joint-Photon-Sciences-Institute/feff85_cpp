// Pure C++ implementations of LAPACK cgetrf_ and cgetrs_.
// These replace the Fortran LAPACK calls in movrlp.cpp and ovp2mt.cpp
// for single-precision complex LU factorization and solve.

#include <complex>
#include <cmath>
#include <algorithm>

using Cf = std::complex<float>;

extern "C" {

// LU factorization with partial pivoting of a complex single-precision matrix.
// Column-major layout, dimensions m x n, leading dimension lda.
// A is overwritten with L (unit lower) and U (upper).
// ipiv contains 1-based pivot indices.
void cgetrf_(const int* m, const int* n, Cf* a,
             const int* lda, int* ipiv, int* info) {
    int rows = *m;
    int cols = *n;
    int ld = *lda;
    *info = 0;

    int minmn = std::min(rows, cols);

    // Column-major access: A(i,j) = a[i + j*ld]
    #define A(i, j) a[(i) + (j) * ld]

    for (int k = 0; k < minmn; ++k) {
        // Find pivot
        int maxidx = k;
        float maxval = std::abs(A(k, k));
        for (int i = k + 1; i < rows; ++i) {
            float val = std::abs(A(i, k));
            if (val > maxval) {
                maxval = val;
                maxidx = i;
            }
        }
        ipiv[k] = maxidx + 1;  // 1-based

        if (maxval == 0.0f) {
            *info = k + 1;
            return;
        }

        // Swap rows k and maxidx
        if (maxidx != k) {
            for (int j = 0; j < cols; ++j) {
                std::swap(A(k, j), A(maxidx, j));
            }
        }

        // Compute multipliers
        Cf pivot = A(k, k);
        for (int i = k + 1; i < rows; ++i) {
            A(i, k) /= pivot;
        }

        // Update trailing submatrix
        for (int j = k + 1; j < cols; ++j) {
            for (int i = k + 1; i < rows; ++i) {
                A(i, j) -= A(i, k) * A(k, j);
            }
        }
    }

    #undef A
}

// Solve A*X = B using the LU factorization from cgetrf_.
// B is overwritten with the solution X.
void cgetrs_(const char* trans, const int* n, const int* nrhs,
             const Cf* a, const int* lda,
             const int* ipiv, Cf* b, const int* ldb,
             int* info) {
    int nn = *n;
    int ld_a = *lda;
    int ld_b = *ldb;
    int nrh = *nrhs;
    *info = 0;

    #define LU(i, j)  a[(i) + (j) * ld_a]
    #define B(i, j)   b[(i) + (j) * ld_b]

    bool notrans = (trans[0] == 'N' || trans[0] == 'n');
    (void)notrans;  // Only 'N' is used in this codebase

    // Apply row permutations (forward): for i=0..n-1, swap row i with ipiv[i]-1
    for (int i = 0; i < nn; ++i) {
        int pi = ipiv[i] - 1;  // convert to 0-based
        if (pi != i) {
            for (int j = 0; j < nrh; ++j) {
                std::swap(B(i, j), B(pi, j));
            }
        }
    }

    // Solve L*Y = P*B (forward substitution, L is unit lower triangular)
    for (int k = 0; k < nn; ++k) {
        for (int i = k + 1; i < nn; ++i) {
            for (int j = 0; j < nrh; ++j) {
                B(i, j) -= LU(i, k) * B(k, j);
            }
        }
    }

    // Solve U*X = Y (back substitution)
    for (int k = nn - 1; k >= 0; --k) {
        for (int j = 0; j < nrh; ++j) {
            B(k, j) /= LU(k, k);
        }
        for (int i = 0; i < k; ++i) {
            for (int j = 0; j < nrh; ++j) {
                B(i, j) -= LU(i, k) * B(k, j);
            }
        }
    }

    #undef LU
    #undef B
}

} // extern "C"
