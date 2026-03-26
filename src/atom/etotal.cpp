// Total energy calculation — converted from src/ATOM/etotal.f
//
// Computes the total atomic energy from:
//   ener[0] = coul  (direct Coulomb Fk integrals)
//   ener[1] = ech   (exchange Coulomb Gk integrals)
//   ener[2] = mag   (magnetic Breit interaction)
//   ener[3] = ret   (retardation Breit interaction)
//
// Note: The Fortran routine receives xnval as a parameter to control
// core-valence separation in the exchange terms. For pure Dirac-Fock
// (idfock=1, the only mode used), xnvalp is all zeros, so the xnval
// checks are never triggered. This implementation matches that behavior.

#include "etotal.hpp"
#include "radial_integrals.hpp"
#include "utility.hpp"
#include "../common/logging.hpp"
#include <feff/constants.hpp>
#include <cmath>
#include <cstdio>

namespace feff::atom {

void etotal(int io, double& eatom, AtomState& state) {
    auto& config = state.config;
    auto& scf    = state.scf;
    auto& breit  = state.breit;
    auto& orb    = state.orb;
    auto& work   = state.work;
    auto& mesh   = state.mesh;
    auto& inel   = state.inelastic;
    auto& ang    = state.angular;
    auto& print_ctl = state.print;

    int norb = scf.norb;

    double ener[4] = {};
    const char* iner[4] = {"coul", "ech.", "mag.", "ret."};

    // Fk integrals (direct Coulomb)
    // Fortran loops i=1,norb; j=1,i  (1-based)
    // C++ loops i=0..norb-1; j=0..i  (0-based)
    int iv = 0;
    for (int i = 0; i < norb; ++i) {
        int l = std::abs(config.kap[i]) - 1;
        for (int j = 0; j <= i; ++j) {
            double a = 1.0;
            if (j == i) a = 2.0;
            int m = std::abs(config.kap[j]) - 1;
            int kmi = 2 * std::min(l, m);
            int k = 0;
            while (k <= kmi) {
                double cer = fdrirk(i, i, j, j, k, orb, work, config, inel, mesh);
                ener[0] += cer * akeato(i, j, k, ang) / a;
                iv++;
                if (iv >= 3) iv = 0;
                k += 2;
            }
        }
    }

    // Gk integrals (exchange Coulomb)
    // For idfock=1: xnvalp is all zeros, so a=1.0 always
    // and the "goto 70" (skip) on xnval(j)>0 never fires.
    iv = 0;
    if (norb > 1) {
        for (int i = 1; i < norb; ++i) {
            double a = 1.0;
            // For idfock != 1: if (xnval[i] > 0.0) a = 0.5;
            for (int j = 0; j < i; ++j) {
                // For idfock != 1: if (xnval[j] > 0.0) continue;  // goto 70
                int l = std::abs(config.kap[i]);
                int m = std::abs(config.kap[j]);
                int k = std::abs(l - m);
                if ((config.kap[i] * config.kap[j]) < 0) k += 1;
                int kmi = l + m - 1;
                while (k <= kmi) {
                    double cer = fdrirk(i, j, i, j, k, orb, work, config, inel, mesh);
                    ener[1] -= cer * bkeato(i, j, k, ang) * a;
                    iv++;
                    if (iv >= 3) iv = 0;
                    k += 2;
                }
            }
        }
    }

    // Breit interaction terms
    inel.nem = 1;

    // Direct Breit integrals
    int ik = 0;
    for (int j = 0; j < norb; ++j) {
        int jj = 2 * std::abs(config.kap[j]) - 1;
        for (int i = 0; i <= j; ++i) {
            int ji = 2 * std::abs(config.kap[i]) - 1;
            int k = 1;
            int kma = std::min(ji, jj);
            while (k <= kma) {
                double cer = fdrirk(j, j, i, i, k, orb, work, config, inel, mesh);
                if (i == j) {
                    bkmrdf(j, j, k, config, breit);
                    ener[2] += (breit.cmag[0] + breit.cmag[1] + breit.cmag[2])
                             * cer * fdmocc(j, j, config) / 2.0;
                }
                ik++;
                if (ik >= 3) ik = 0;
                k += 2;
            }
        }
    }

    // Exchange Breit integrals
    if (norb > 1) {
        for (int j = 1; j < norb; ++j) {
            int lj = std::abs(config.kap[j]);
            int na = -1;
            if (config.kap[j] <= 0) {
                na = 1;
                lj = lj - 1;
            }
            for (int l = 0; l < j; ++l) {
                int ll = std::abs(config.kap[l]);
                int nb = -1;
                if (config.kap[l] <= 0) {
                    nb = 1;
                    ll = ll - 1;
                }
                double b = fdmocc(j, l, config);
                int nm1  = std::abs(lj + na - ll);
                int nmp1 = ll + lj + nb;
                int nmm1 = ll + lj + na;
                int np1  = std::abs(ll + nb - lj);
                int k = std::min(nm1, np1);
                int kma = std::max(nmp1, nmm1);
                if ((k + ll + lj) % 2 == 0) k += 1;
                int nb_sum = std::abs(config.kap[j]) + std::abs(config.kap[l]);

                while (k <= kma) {
                    bkmrdf(j, l, k, config, breit);
                    double cer[3] = {};

                    // Fortran: if (nb.le.k.and.kap(l).lt.0.and.kap(j).gt.0) go to 161
                    if (!(nb_sum <= k && config.kap[l] < 0 && config.kap[j] > 0)) {
                        cer[0] = fdrirk(l, j, l, j, k, orb, work, config, inel, mesh);
                        // Fortran passes 0,0 meaning "reuse workspace yk"
                        // In 0-based C++, we pass -1,-1 as sentinel
                        cer[1] = fdrirk(-1, -1, j, l, k, orb, work, config, inel, mesh);
                    }

                    // Fortran: if (nb.le.k.and.kap(l).gt.0.and.kap(j).lt.0) go to 171
                    if (!(nb_sum <= k && config.kap[l] > 0 && config.kap[j] < 0)) {
                        cer[2] = fdrirk(j, l, j, l, k, orb, work, config, inel, mesh);
                        if (cer[1] == 0.0) {
                            cer[1] = fdrirk(-1, -1, l, j, k, orb, work, config, inel, mesh);
                        }
                    }

                    for (int ii = 0; ii < 3; ++ii) {
                        ener[2] += breit.cmag[ii] * cer[ii] * b;
                        ener[3] += breit.cret[ii] * cer[ii] * b;
                    }
                    k += 2;
                }
            }
        }
    }

    // Total energy: E = -(Coul + Exch) + Mag + Ret + sum(en * xnel)
    eatom = -(ener[0] + ener[1]) + ener[2] + ener[3];
    for (int j = 0; j < norb; ++j) {
        eatom += config.en[j] * config.xnel[j];
    }

    // Logging output
    if (print_ctl.iprint >= 5) {
        char buf[128];
        std::snprintf(buf, sizeof(buf), "etot %18.7e", eatom * feff::hart);
        feff::common::logger().wlog(buf);
        for (int i = 0; i < 4; ++i) {
            std::snprintf(buf, sizeof(buf), "%s %18.7e", iner[i], ener[i] * feff::hart);
            feff::common::logger().wlog(buf);
        }
    }
}

} // namespace feff::atom
