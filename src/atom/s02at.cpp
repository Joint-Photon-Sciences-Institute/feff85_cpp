// S02 overlap integral calculation — converted from src/ATOM/s02at.f
//
// Calculates the S0^2 shakeup amplitude from orbital overlap integrals.
// Uses determinant formulation for multi-orbital overlap between
// initial and final state configurations.

#include "s02at.hpp"
#include "../math/determinant.hpp"
#include <cmath>
#include <cstring>

namespace feff::atom {

void elimin(double* d1, int n, double* d2) {
    // d1 and d2 are 7x7 matrices (row-major)
    // n is 1-based row/column to replace with identity row/column
    // In Fortran: d1(7,7), d2(7,7) with column-major
    // d1(i,j) → d1[(j-1)*7 + (i-1)] in column-major
    // We use row-major: d1[i*7 + j] where i,j are 0-based
    //
    // But to match the Fortran determ() which expects a specific layout,
    // we store in column-major compatible order: d1[j*7 + i] for (i,j) 1-based
    // Actually determ takes (array, nord, nrows) where nrows is leading dim.
    // The Fortran stores column-major. Our determ.cpp should match.
    //
    // For consistency, treat matrices as d[i][j] with i=row, j=col, 0-based.
    // Access: d[i*7 + j]
    // n is 1-based index in Fortran; convert to 0-based
    int n0 = n - 1;

    for (int i = 0; i < 7; ++i) {
        for (int j = 0; j < 7; ++j) {
            if (i != n0) {
                if (j != n0) {
                    d2[i * 7 + j] = d1[i * 7 + j];
                } else {
                    d2[i * 7 + j] = 0.0;
                }
            } else {
                if (j != n0) {
                    d2[i * 7 + j] = 0.0;
                } else {
                    d2[i * 7 + j] = 1.0;
                }
            }
        }
    }
}

void s02at(int ihole, int norb, const int nk[30], const double xnel[30],
           const double ovpint[][30], double& dval)
{
    // ihole is 0-based orbital index (iholep from Fortran, already converted)
    // nk[i] = kappa quantum number for orbital i (0-based)
    // ovpint[i][j] = overlap integral between initial orbital i and final orbital j

    dval = 1.0;

    // Loop over possible kappa values: -4 to 3
    for (int kap = -4; kap <= 3; ++kap) {
        // Initialize 7x7 matrices
        double m1[7 * 7] = {};
        double m2[7 * 7] = {};
        int iorb[7] = {};

        for (int i = 0; i < 7; ++i) {
            m1[i * 7 + i] = 1.0;  // identity
            m2[i * 7 + i] = 1.0;
        }

        int morb = 0;
        int nhole = 0;

        // Construct the overlap matrix for orbitals with this kappa
        for (int i = 0; i < norb; ++i) {
            if (nk[i] == kap) {
                iorb[morb] = i;
                // Fill row/column morb with overlap integrals
                for (int j = 0; j <= morb; ++j) {
                    m1[j * 7 + morb] = ovpint[iorb[j]][iorb[morb]];
                }
                for (int j = 0; j < morb; ++j) {
                    m1[morb * 7 + j] = m1[j * 7 + morb];
                }

                if (ihole == i) nhole = morb + 1;  // 1-based nhole as in Fortran
                morb++;
            }
        }

        if (morb != 0) {
            // Make working copies for determinant (it modifies in place)
            double m1_copy[7 * 7];
            double m1_copy2[7 * 7];

            std::memcpy(m1_copy, m1, sizeof(m1));
            double dum1 = feff::math::determ(m1_copy, morb, 7);
            dum1 = dum1 * dum1;

            std::memcpy(m1_copy2, m1, sizeof(m1));
            double dum3 = feff::math::determ(m1_copy2, morb - 1, 7);
            dum3 = dum3 * dum3;

            // xnel uses iorb[morb-1] (last orbital added)
            double xn = xnel[iorb[morb - 1]];
            int nmax_kap = 2 * std::abs(kap);
            double xnh = nmax_kap - xn;

            if (nhole == 0) {
                dval *= std::pow(dum1, xn) * std::pow(dum3, xnh);
            } else if (nhole == morb) {
                // nhole is 1-based, morb is count
                dval *= std::pow(dum1, xn - 1) * std::pow(dum3, xnh + 1);
            } else {
                // nhole is 1-based index into the morb-sized submatrix
                double m1_copy3[7 * 7];
                double m2_copy[7 * 7];
                double m2_copy2[7 * 7];

                std::memcpy(m1_copy3, m1, sizeof(m1));
                elimin(m1_copy3, nhole, m2);

                std::memcpy(m2_copy, m2, sizeof(m2));
                double dum2 = feff::math::determ(m2_copy, morb, 7);
                dum2 = dum2 * dum2;

                std::memcpy(m2_copy2, m2, sizeof(m2));
                double dum4 = feff::math::determ(m2_copy2, morb - 1, 7);
                dum4 = dum4 * dum4;

                double dum5 = (dum4 * dum1 * xnh + dum2 * dum3 * xn) / nmax_kap;
                dval *= dum5 * std::pow(dum1, xn - 1) * std::pow(dum3, xnh - 1);
            }
        }
    }
}

} // namespace feff::atom
