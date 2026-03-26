// Alternative energy mesh with user-defined grid support.
// Converted from src/XSPH/phmesh2.f
//
// Contains helper subroutines: MkVGrid84, MkExpMesh, ExafsGrid84,
// XanesGrid84, FPrimeGrid84, ReverseGrid, MkEMesh, MkKMesh, SortE

#include "phmesh2.hpp"
#include "rdgrid.hpp"
#include <feff/constants.hpp>
#include <feff/dimensions.hpp>
#include "../common/qsort.hpp"
#include <cmath>
#include <algorithm>
#include <fstream>

namespace feff::common { double getxk(double e); }

namespace feff::xsph {

// --- Helper routines (local to this translation unit) ---

static void mk_emesh(FeffComplex em[], int iStart, double emin, double emax,
                     double estep, int& npts, int nex_dim) {
    npts = static_cast<int>(std::round((emax - emin) / estep));
    if (npts <= 0) { npts = 0; return; }
    for (int i = 0; i <= npts; i++) {
        if (iStart + i < nex_dim) {
            em[iStart + i] = emin + estep * i;
        }
    }
}

static void mk_kmesh(FeffComplex em[], int iStart, double xkmin, double xkmax,
                     double deltak, int& npts, int nex_dim) {
    npts = static_cast<int>(std::round((xkmax - xkmin) / deltak));
    if (npts <= 0) { npts = 0; return; }
    int isgn = (xkmin < 0.0) ? -1 : 1;
    for (int i = 0; i <= npts; i++) {
        if (iStart + i < nex_dim) {
            double k = xkmin + deltak * i;
            em[iStart + i] = isgn * k * k / 2.0;
        }
    }
}

static void mk_expmesh(FeffComplex em[], int iStart, double emin, double emax,
                       double del, int& npts, int nex_dim) {
    npts = static_cast<int>(std::round(std::log(emax / emin) / del));
    for (int i = 0; i <= npts; i++) {
        if (iStart + i < nex_dim) {
            em[iStart + i] = emin * std::exp(del * i);
        }
    }
}

static void exafs_grid84(FeffComplex em[], double xkmax, int& ne, int nex_dim) {
    int npts;
    // 20 pts (0 <= k <= 1.9, delk=0.1/A)
    double deltak = bohr / 10.0;
    mk_kmesh(em, ne, 0.0, bohr * 1.9 * 1.01, deltak, npts, nex_dim);
    // 20 pts (2 <= k <= 5.8, delk=0.2/A)
    ne += npts + 1;
    deltak = bohr / 5.0;
    mk_kmesh(em, ne, bohr * 2.0, bohr * 5.8 * 1.01, deltak, npts, nex_dim);
    // 9 pts (6 <= k <= 10, delk=0.5/A)
    ne += npts + 1;
    mk_kmesh(em, ne, bohr * 6.0, bohr * 10.0 * 1.01, bohr * 0.5, npts, nex_dim);
    // rest with deltak=1.0/A
    ne += npts + 1;
    deltak = bohr;
    double xkmin = std::sqrt(2.0 * std::real(em[ne - 1])) + deltak;
    npts = std::min(100 - ne, static_cast<int>(std::round((xkmax - xkmin) / deltak)) + 1);
    double xkmax2 = xkmin + npts * deltak * 1.01;
    mk_kmesh(em, ne, xkmin, xkmax2, deltak, npts, nex_dim);
    ne += npts;
}

static void xanes_grid84(FeffComplex em[], double xkmax, double xkstep,
                         double estep, int& ne, int& ik0, int nex_dim) {
    int nemax = 10;
    double dk = 2.0 * xkstep;
    int n1 = static_cast<int>(estep / 2.0 / (dk * dk));
    int n2 = static_cast<int>(std::sqrt(n1 * 2.0 * estep) / dk);
    if ((dk * (n2 + 1)) * (dk * (n2 + 1)) > (n1 + 1) * 2.0 * estep) n1++;
    n1 = std::min(n1, nemax);
    int nk = nemax - n1;

    // Fill k grid below
    int npts;
    double xkmin = -dk * (n2 + nk);
    double xkmax2 = -dk * (n2 + 1);
    ne = 0;
    mk_kmesh(em, ne, xkmin, xkmax2, dk, nk, nex_dim);

    // Fill e grid
    ne += nk + 1;
    double emin = -estep * n1;
    mk_emesh(em, ne, emin, 0.0, estep, npts, nex_dim);
    ne += npts + 1;
    ik0 = ne;

    // Above Fermi level
    nemax = 90;
    n1 = static_cast<int>(estep / 2.0 / (xkstep * xkstep));
    n2 = static_cast<int>(std::sqrt(n1 * 2.0 * estep) / xkstep);
    n1++;
    if ((xkstep * (n2 + 1)) * (xkstep * (n2 + 1)) > n1 * 2.0 * estep) n1++;
    n1 = std::min(n1, nemax);
    nk = nemax - n1;

    emin = estep;
    double emax = (n1 - 1) * estep;
    if (std::sqrt(2.0 * emax) > xkmax) {
        emax = xkmax * xkmax / 2.0;
        nk = 0;
    }
    mk_emesh(em, ne, emin, emax, estep, npts, nex_dim);

    ne += npts + 1;
    xkmin = xkstep * (n2 + 1);
    xkmax2 = xkstep * (n2 + nk);
    if (xkmax2 > xkmax) xkmax2 = xkmax;
    mk_kmesh(em, ne, xkmin, xkmax2, xkstep, npts, nex_dim);
    ne += npts;
}

static void fprime_grid84(FeffComplex em[], double emin_in, double emax_in,
                          double estep, double emu, double edge,
                          int& ne, int& ne1, int& ne3, int nex_dim) {
    int nemax = 100;
    double emin = emin_in / bohr / hart - emu;
    double emax = emax_in / bohr / hart - emu;

    em[0] = emin;
    ne = 1;
    if (emin < emax) {
        if (estep <= 0.0) estep = (emax - emin) / (nemax - 1);
        ne = std::min(nemax, static_cast<int>(std::round((emax - emin) / estep)));
        for (int i = 0; i < ne; i++) {
            em[i] = emin + (i + 1) * estep;
        }
    }
    ne1 = ne;

    // KK-Transform grid
    nemax = std::min(nex_dim - ne, 100);
    double del = 3.0 / hart;
    double elimit = std::max(1.0e3 / hart, std::min(20.0 * emu, 2.0e5 / hart));
    elimit = elimit - emu;

    ne3 = nemax;
    em[ne1] = edge;
    for (int i = 1; i < ne3; i++) {
        double del2 = 0.0;
        if (std::real(em[ne1 + i - 1]) > 0.0) {
            del2 = std::real(em[ne1 + i - 1]) *
                (std::exp(std::log(elimit / std::real(em[ne1 + i - 1])) / (ne3 - i)) - 1.0);
        }
        em[ne1 + i] = em[ne1 + i - 1] + std::max(del, del2);
    }
    ne = ne1 + ne3;
}

static void mk_vgrid84(FeffComplex em[], int& ne, double xloss, int nex_dim) {
    double estep0 = 0.01 / hart;
    em[ne] = coni * estep0 / 2.0;
    em[ne + 1] = coni * estep0;
    ne += 2;

    double del = 0.4;
    int n1 = static_cast<int>(std::round(std::log(xloss / estep0) / del - 0.5));
    if (n1 <= 0) n1 = 1;
    double expdel = std::exp(del);
    double emin = 2.0 * xloss / (1.0 + expdel) / std::pow(expdel, n1);
    if (emin <= estep0) emin = emin * expdel;
    double emax = std::min(50.0 / hart, 20.0 * xloss);
    int npts;
    mk_expmesh(em, ne, emin, emax, del, npts, nex_dim);
    for (int i = 0; i <= npts; i++) {
        em[ne + i] = coni * em[ne + i];
    }
    ne += npts;
}

static void reverse_grid(FeffComplex em[], int ne, double zero_point) {
    for (int i = 0; i < ne; i++) {
        em[i] = zero_point - em[i];
    }
    int np = ne / 2;
    for (int i = 0; i < np; i++) {
        FeffComplex tmp = em[i];
        em[i] = em[ne - 1 - i];
        em[ne - 1 - i] = tmp;
    }
}

static void sort_e(FeffComplex em[], int& ne, int& ik0, int nex_dim) {
    double tol = 0.001;
    std::vector<double> realE(ne);
    for (int i = 0; i < ne; i++) realE[i] = std::real(em[i]);
    auto order = feff::common::argsort(realE);

    // Replace em with sorted values and remove degeneracy
    std::vector<FeffComplex> sorted(ne);
    int nue = 1;
    if (std::abs(realE[order[0]]) < tol) {
        sorted[0] = 0.0;
        ik0 = 0;
    } else {
        sorted[0] = realE[order[0]];
    }
    for (int i = 1; i < ne; i++) {
        if (std::abs(realE[order[i]] - std::real(sorted[nue - 1])) > tol) {
            sorted[nue] = realE[order[i]];
            nue++;
        }
    }
    ne = nue;
    for (int i = 0; i < ne; i++) em[i] = sorted[i];

    // Set ik0
    ik0 = 0;
    double e0 = std::abs(std::real(em[0]));
    for (int i = 0; i < ne; i++) {
        if (std::abs(std::real(em[i])) < e0) {
            e0 = std::abs(std::real(em[i]));
            ik0 = i;
        }
    }
    em[ik0] = 0.0;
}

// --- Main function ---

void phmesh2(int iprint, int ispec, double edge, double emu,
             double vi0, double gamach,
             double xkmax, double xkstep, double vixan,
             int& ne, int& ne1, FeffComplex em[], int& ik0, int& ne3,
             int iGrid) {

    double xloss = std::max(gamach / 2.0 + vi0, 0.02 / hart);
    double xim = (vixan > 0.0001) ? vixan : xloss / 2.0;
    ik0 = 0;

    if (iGrid == 0) {
        // Use FEFF84 grids
        if (ispec == 0) {
            ne = 1;
            exafs_grid84(em, xkmax, ne, nex);
            ne1 = ne;
            ik0 = 0;
        } else if (std::abs(ispec) > 0 && std::abs(ispec) < 4) {
            xanes_grid84(em, xkmax, xkstep, xim, ne, ik0, nex);
            ne1 = ne;
        } else if (ispec == 4) {
            fprime_grid84(em, xkmax, xkstep, vixan, emu, edge, ne, ne1, ne3, nex);
        }
        if (ispec < 0) {
            ne = 11;
            exafs_grid84(em, xkmax, ne, nex);
            ne1 = ne;
        }
    } else {
        // User defined grids
        int nemax = nex - 50;
        ne = 0;
        constexpr int nGridMax = 10;
        int nGrid;
        int iGridType[nGridMax];
        double gridMin[nGridMax], gridMax[nGridMax], gridStep[nGridMax];
        rdgrid(em, ne, nGrid, iGridType, gridMin, gridMax, gridStep, nGridMax, nemax);

        for (int g = 0; g < nGrid; g++) {
            int npts;
            if (iGridType[g] == 1) {
                ne++;
                mk_emesh(em, ne, gridMin[g], gridMax[g], gridStep[g], npts, nex);
                ne = std::min(ne + npts, nemax);
            } else if (iGridType[g] == 2) {
                ne++;
                mk_kmesh(em, ne, gridMin[g], gridMax[g], gridStep[g], npts, nex);
                ne = std::min(ne + npts, nemax);
            } else if (iGridType[g] == 3) {
                ne++;
                mk_expmesh(em, ne, gridMin[g], gridMax[g], gridStep[g], npts, nex);
                ne = std::min(ne + npts, nemax);
            }
        }
        // Add a point at E=0
        if (ne + 1 < nex) {
            em[ne] = 0.0;
            ne++;
        } else {
            em[ne - 1] = 0.0;
        }
        sort_e(em, ne, ik0, nex);
        ne1 = ne;
    }

    // If XES, flip grid about 0.0
    if (std::abs(ispec) == 2) reverse_grid(em, ne, 0.0);

    // Shift horizontal grid by edge + i*xloss
    if (ispec != 4) {
        for (int i = 0; i < ne; i++) {
            em[i] = em[i] + edge + coni * xloss;
        }
    }

    // If not fprime, make vertical grid
    if (ispec != 4) {
        ne++;
        mk_vgrid84(em, ne, xloss, nex);
        // Shift vertical grid by edge
        for (int i = ne1; i < ne; i++) {
            em[i] = em[i] + edge;
        }
    }

    if (std::abs(ispec) == 3) {
        // DANES: add more points
        ne3 = std::min(nex, 150) - ne;
        double emin = std::real(2.0 * em[ne1 - 1] - em[ne1 - 2]);
        double emax = 7.0e4;
        double del = std::log(emax / emin) / (ne3 - 1);
        ne++;
        int npts;
        mk_expmesh(em, ne, emin, emax, del, npts, nex);
        for (int i = 0; i < ne3; i++) {
            em[ne + i] = em[ne + i] + coni * 1.0e-8;
        }
        ne = ne + ne3;
    }

    if (iprint >= 3) {
        std::ofstream fout("emesh.dat");
        fout << "edge, bohr, edge*hart " << edge << " " << bohr << " " << edge * hart << "\n";
        fout << "ispec, ik0 " << ispec << " " << ik0 << "\n";
        fout << "ie, em(ie)*hart, xk(ie)\n";
        for (int ie = 0; ie < ne; ie++) {
            fout << ie << " " << std::real(em[ie]) * hart << " "
                 << feff::common::getxk(std::real(em[ie]) - edge) / bohr << "\n";
        }
    }
}

} // namespace feff::xsph
