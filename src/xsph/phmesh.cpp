// Energy mesh generation for phase shift calculations.
// Converted from src/XSPH/phmesh.f

#include "phmesh.hpp"
#include <feff/constants.hpp>
#include <feff/dimensions.hpp>
#include <cmath>
#include <algorithm>
#include <fstream>

namespace feff::common { double getxk(double e); }

namespace feff::xsph {

void phmesh(int iprint, int ispec, double edge, double emu,
            double vi0, double gamach,
            double xkmax, double xkstep, double vixan,
            int& ne, int& ne1, FeffComplex em[], int& ik0, int& ne3) {

    double xloss = gamach / 2.0 + vi0;
    if (xloss < 0.0) xloss = 0.0;
    double xvert = std::max(xloss, 0.02 / hart);
    xloss = xvert;
    double aa = 0.5;
    ne3 = 0;
    double xim = xloss * aa;
    if (vixan > 0.0001) xim = vixan;
    ik0 = 0;

    int nmin = 0;

    if (ispec <= 3) {
        // Make energy mesh for XANES with FMS
        // 10 points below Fermi level
        int nemax = 10;
        double dk = 2.0 * xkstep;
        int n1 = static_cast<int>(xim / 2.0 / (dk * dk));
        int n2 = static_cast<int>(std::sqrt(n1 * 2.0 * xim) / dk);
        if ((dk * (n2 + 1)) * (dk * (n2 + 1)) > (n1 + 1) * 2.0 * xim) n1++;
        n1 = std::min(n1, nemax);
        for (int i = 1; i <= n1; i++) {
            em[nemax - i] = -xim * i + edge + coni * xloss;
        }
        int nb = nemax - n1;
        for (int i = 1; i <= nb; i++) {
            em[nb - i] = -(dk * (n2 + i)) * (dk * (n2 + i)) / 2.0 + edge + coni * xloss;
        }
        nmin = nemax;
        ik0 = nemax; // 0-based index where k=0
    }

    if (ispec > 0 && ispec <= 3) {
        // 90 points above Fermi level
        int nemax = 100 - 10;
        int n1 = static_cast<int>(xim / 2.0 / (xkstep * xkstep));
        int n2 = static_cast<int>(std::sqrt(n1 * 2.0 * xim) / xkstep);
        n1++;
        if ((xkstep * (n2 + 1)) * (xkstep * (n2 + 1)) > n1 * 2.0 * xim) n1++;
        n1 = std::min(n1, nemax);
        int nb;
        if (ispec != 2) {
            nb = static_cast<int>(xkmax * xkmax / xim / 2.0) + 1;
        } else {
            nb = static_cast<int>(std::abs(edge - xkmax / bohr / hart) / xim) + 1;
        }
        if (nb <= n1) n1 = nb;
        for (int i = 1; i <= n1; i++) {
            em[nmin + i - 1] = xim * (i - 1);
        }
        if (ispec != 2) {
            nb = static_cast<int>(xkmax / xkstep) - n2;
        } else {
            nb = static_cast<int>(std::sqrt(std::abs(2.0 * (edge - xkmax / bohr / hart))) / xkstep) - n2;
        }
        nb = std::min(nb, nemax - n1);
        nb = std::max(nb, 0);
        for (int i = 1; i <= nb; i++) {
            em[nmin + n1 + i - 1] = (xkstep * (n2 + i)) * (xkstep * (n2 + i)) / 2.0;
        }
        ne1 = nmin + n1 + nb;
        for (int i = ik0; i < ne1; i++) {
            em[i] = em[i] + edge + coni * xloss;
        }

    } else if (ispec == 4) {
        // Grid for atomic f' calculation
        int nemax = 100;
        double emin = xkmax / bohr / hart;
        double emax = xkstep / bohr / hart;
        ne = 1;
        emin = emin - emu + edge;
        emax = emax - emu + edge;
        em[0] = emin;
        if (emin < emax) {
            if (vixan <= 0.0) vixan = (emax - emin) / (nemax - 1);
            while (ne < nemax && std::real(em[ne - 1]) < emax) {
                em[ne] = em[ne - 1] + vixan;
                ne++;
            }
        }
        ne1 = ne;
        nemax = nex - ne;
        if (nemax > 100) nemax = 100;
        double de = 3.0 / hart;
        double elimit = std::min(2.0e5 / hart, 20.0 * emu);
        elimit = std::max(elimit, 1.0e3 / hart);
        elimit = elimit - emu;

        int ne2 = 0;
        ne3 = nemax;
        ne = ne1 + ne2 + ne3;
        em[ne1] = edge;
        for (int i = 1; i < ne3; i++) {
            double dep = 0.0;
            if (std::real(em[ne1 + i - 1]) > 0.0) {
                dep = std::real(em[ne1 + i - 1]) *
                    (std::exp(std::log(elimit / std::real(em[ne1 + i - 1])) / (ne3 - i)) - 1.0);
            }
            if (dep < de) dep = de;
            em[ne1 + i] = em[ne1 + i - 1] + dep;
        }
    } else {
        // Energy mesh for EXAFS or XANES without FMS
        ne = 0;
        if (ispec < 0) ne = 10;
        double delk = bohr / 10.0;
        for (int i = 0; i < 20; i++) {
            double tempk = i * delk;
            em[ne] = tempk * tempk / 2.0 + edge + coni * xloss;
            if (i == 0) ik0 = ne;
            ne++;
        }
        delk = bohr / 5.0;
        for (int i = 0; i < 20; i++) {
            double tempk = 2.0 * bohr + i * delk;
            em[ne] = tempk * tempk / 2.0 + edge + coni * xloss;
            ne++;
        }
        delk = bohr / 2.0;
        for (int i = 0; i < 9; i++) {
            double tempk = 6.0 * bohr + i * delk;
            em[ne] = tempk * tempk / 2.0 + edge + coni * xloss;
            ne++;
        }
        delk = bohr;
        double tempk = 0.0;
        for (int i = 0; i < 10; i++) {
            tempk = 11.0 * bohr + i * delk;
            em[ne] = tempk * tempk / 2.0 + edge + coni * xloss;
            ne++;
        }
        while (tempk < xkmax) {
            tempk = tempk + delk;
            em[ne] = tempk * tempk / 2.0 + edge + coni * xloss;
            ne++;
        }
        ne = std::min(ne, 100);
        ne1 = ne;
    }

    if (ispec <= 3) {
        // Make the vertical grid in energy plane
        double tempk = 0.005 / hart;
        em[ne1] = edge + coni * tempk;
        tempk = tempk * 2.0;
        em[ne1 + 1] = edge + coni * tempk;
        double delk = 0.4;
        int n1 = static_cast<int>(std::round(std::log(xloss / tempk) / delk - 0.5));
        if (n1 <= 0) n1 = 1;
        double bb = std::exp(delk);
        aa = 2.0 * xloss / (1.0 + bb);
        aa = aa / std::pow(bb, n1);
        if (aa <= tempk) aa = aa * bb;
        double ee = std::min(50.0 / hart, 20.0 * xloss);
        n1 = static_cast<int>(std::round(std::log(ee / aa) / delk));
        for (int i = 0; i <= n1; i++) {
            em[ne1 + 2 + i] = edge + coni * aa * std::exp(delk * i);
        }
        ne = ne1 + n1 + 3;

        // For DANES: additional points
        if (std::abs(ispec) == 3) {
            ne3 = std::min(nex, 150) - ne;
            em[ne] = std::real(2.0 * em[ne1 - 1] - em[ne1 - 2]);
            double dk2 = std::log(7.0e4 / std::real(em[ne])) / (ne3 - 1);
            dk2 = std::exp(dk2);
            for (int i = 1; i < ne3; i++) {
                em[ne + i] = em[ne + i - 1] * dk2;
            }
            for (int i = 0; i < ne3; i++) {
                em[ne + i] = em[ne + i] + coni * 1.0e-8;
            }
            ne = ne + ne3;
        }
    }

    // Reverse order for horizontal grid for XES
    if (ispec == 2) {
        for (int ie = 0; ie < ne1; ie++) {
            em[ie] = 2.0 * (edge + coni * xloss) - em[ie];
        }
        int np = ne1 / 2;
        for (int ie = 0; ie < np; ie++) {
            int ip = ne1 - 1 - ie;
            FeffComplex tempc = em[ie];
            em[ie] = em[ip];
            em[ip] = tempc;
        }
        ik0 = ne1 - 1 - ik0;
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
