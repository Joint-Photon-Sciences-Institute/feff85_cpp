// Core-valence separation.
// Converted from src/POT/corval.f

#include "corval.hpp"
#include "rholie.hpp"
#include "../common/grid_interpolation.hpp"
#include "../common/logging.hpp"
#include <feff/constants.hpp>
#include <cmath>
#include <complex>
#include <sstream>

namespace feff::pot {

void corval(bool verbse,
            double& ecv, double* xnvmu, const double* eorb, const int* norb,
            double* xnval, const int* kappa, double rgrd,
            int nohole, int nph, double* edens, double* edenvl,
            const double* vtot, const double* vvalgs,
            const double* rmt, const double* rnrm, int ixc,
            double rhoint, double vint, int jumprm,
            double x0, double* ri, double dx, const double* xion,
            int iunf, const int* iz,
            const double* adgc, const double* adpc,
            double* dgc, double* dpc,
            int ihole, const int* lmaxsc)
{
    using FeffComplex = std::complex<double>;
    constexpr int s251 = 251;

    // Local arrays — use static for large arrays to avoid stack overflow.
    // corval calls rholie which also has large arrays; combined they exceed stack.
    static double dmagx[nrptx], dmag0[s251];
    static double vtotph[nrptx], vvalph[nrptx];
    static double dum[nrptx];
    static double ri05[s251];
    static double dgcn[nrptx * 30], dpcn[nrptx * 30];
    FeffComplex xrhoce[lx + 1], xrhole[lx + 1];
    FeffComplex yrhoce[s251], yrhole[s251 * (lx + 1)];
    FeffComplex ph[lx + 1];

    // Sorted level arrays (max 32 = (lx+1)*(nphx+1))
    constexpr int max_levels = 32;
    double en[max_levels] = {};
    int ll[max_levels] = {};
    int ip[max_levels] = {};
    int icv[max_levels] = {};

    // Per-potential, per-l arrays
    double eldos[(lx + 1) * (nphx + 1)] = {};
    int iiorb[(lx + 1) * (nphx + 1)] = {};
    int ival[(lx + 1) * (nphx + 1)] = {};
    int ifound[lx + 1];

    if (verbse) {
        feff::common::logger().wlog("              Core-valence separation ");
    }

    // Initialize
    for (int i = 0; i < s251; ++i) {
        dmag0[i] = 0.0;
        ri05[i] = std::exp(-8.8 + 0.05 * i);
    }
    for (int iph = 0; iph <= nphx; ++iph) {
        for (int il = 0; il <= lx; ++il) {
            int idx = iph * (lx + 1) + il;
            eldos[idx] = 0.0;
            iiorb[idx] = 0;
            ival[idx] = 0;
        }
    }

    double tol = 5.0 / hart;
    if (vint - ecv < tol) ecv = vint - tol;
    double elow = -70.0 / hart;
    double ehigh = -20.0 / hart;
    FeffComplex eimag = coni * 1.5 / hart;
    // Make energy step about 0.5 eV
    int ne = 1 + static_cast<int>(std::round((ehigh - elow) * 2.0 * hart));
    double de = (ehigh - elow) / (ne - 1);

    // Find problematic energies for core-valence separation
    for (int iph = 0; iph <= nph; ++iph) {
        for (int iorb = 0; iorb < norb[iph]; ++iorb) {
            double eorb_val = eorb[iph * 30 + iorb];
            if (eorb_val < ehigh - tol && eorb_val > elow) {
                int lll = -kappa[iph * 30 + iorb] - 1;
                if (lll < 0) lll = kappa[iph * 30 + iorb];
                // Skip f-electrons for Hf,Lu,Ta unless UNFREEZEF
                if ((iz[iph] >= 71 && iz[iph] <= 73) && lll == 3) continue;
                if (iunf == 0 && lll == 3) continue;

                int idx = iph * (lx + 1) + lll;
                eldos[idx] = eorb_val;
                ival[idx] = 1;
                if (xnval[iph * 30 + iorb] < 0.1) ival[idx] = -1;
                iiorb[idx] = iorb;
            }
        }
    }

    // Search for suspicious maxima in LDOS for each potential
    for (int iph = 0; iph <= nph; ++iph) {
        double vjump = 0.0;
        feff::common::fixvar(rmt[iph], &edens[iph * s251], &vtot[iph * s251], dmag0,
                             vint, rhoint, dx, rgrd, jumprm,
                             vjump, ri, vtotph, dum, dmagx);
        if ((ixc % 10) >= 5) {
            int jumprm_local = (jumprm > 0) ? 2 : jumprm;
            feff::common::fixvar(rmt[iph], &edenvl[iph * s251], &vvalgs[iph * s251],
                                 dmag0, vint, rhoint, dx, rgrd, jumprm_local,
                                 vjump, ri, vvalph, dum, dmagx);
        }
        feff::common::fixdsx(iph, dx, rgrd, dgc, dpc, dgcn, dpcn);

        int jri = static_cast<int>((std::log(rmt[iph]) + x0) / rgrd) + 2;
        int jri1 = jri + 1;
        // jri1 is 1-based Fortran index => 0-based: jri1-1
        FeffComplex eref = vtotph[jri1 - 1];
        for (int i = 0; i < jri1; ++i) {
            vtotph[i] = vtotph[i] - std::real(eref);
        }
        if (ixc >= 5) {
            for (int i = 0; i < jri1; ++i) {
                vvalph[i] = vvalph[i] - std::real(eref);
            }
        } else {
            for (int i = 0; i < jri1; ++i) {
                vvalph[i] = vtotph[i];
            }
        }

        int itmp = 0;
        if (iph == 0 && nohole < 0) itmp = ihole;

        double xx_val = std::imag(eimag);
        int nfound = 0;
        double xpeak[lx + 1];
        double xp[lx + 1];
        for (int il = 0; il <= lx; ++il) {
            xpeak[il] = (2.0 * il + 1.0) / (6.0 * xx_val * pi);
            xp[il] = 0.0;
            ifound[il] = 1;
            if (ival[iph * (lx + 1) + il] != 0) ifound[il] = 0;
            nfound += ifound[il];
        }
        if (nfound == lx + 1) continue;

        // Start the search for suspicious maxima in LDOS
        int nr05;
        int ie = 0;
        while (ie < ne) {
            ie++;
            FeffComplex emg = elow + de * (ie - 1) + eimag;

            rholie(ri05, nr05, rgrd, x0, ri, emg,
                   ixc, rmt[iph], rnrm[iph],
                   vtotph, vvalph, &xnval[iph * 30], dgcn, dpcn, eref,
                   &adgc[iph * 10 * 30], &adpc[iph * 10 * 30],
                   xrhole, xrhoce, yrhole, yrhoce, ph,
                   iz[iph], xion[iph], iunf, itmp, lmaxsc[iph]);

            // Find suspicious peaks in LDOS
            nfound = 0;
            for (int il = 0; il <= lx; ++il) {
                if (ival[iph * (lx + 1) + il] != 0 && ifound[il] == 0) {
                    double xx_imag = std::imag(xrhoce[il]);
                    if ((ie == ne || xx_imag < xp[il]) && xp[il] > xpeak[il]) {
                        ifound[il] = 1;
                        eldos[iph * (lx + 1) + il] = elow + de * (ie - 2);
                    } else {
                        xp[il] = xx_imag;
                    }
                }
                nfound += ifound[il];
            }
            if (nfound >= lx + 1 || ie >= ne) break;
        }

        if (nfound < lx + 1) {
            feff::common::logger().wlog("WARNING: fatal error in subroutine corval. Try");
            feff::common::logger().wlog("  to reduce ca1 in SCF card. If does not help,");
            feff::common::logger().wlog("SEND bug report to AUTHORS");
            throw std::runtime_error("CORVAL-1");
        }
    }

    // Arrange suspicious levels in order
    ne = 0;
    for (int iph = 0; iph <= nph; ++iph) {
        for (int il = 0; il <= lx; ++il) {
            int idx = iph * (lx + 1) + il;
            if (eldos[idx] < 0.0) {
                ne++;
                // Find insertion position
                int inew = ne;
                for (int ie = 1; ie < ne; ++ie) {
                    if (en[ie - 1] > eldos[idx] && inew == ne) inew = ie;
                }
                // Shift elements right
                for (int ie = ne - 1; ie >= inew; --ie) {
                    en[ie] = en[ie - 1];
                    icv[ie] = icv[ie - 1];
                    ll[ie] = ll[ie - 1];
                    ip[ie] = ip[ie - 1];
                }
                en[inew - 1] = eldos[idx];
                icv[inew - 1] = ival[idx];
                ll[inew - 1] = il;
                ip[inew - 1] = iph;
            }
        }
    }

    // If no suspicious points, exit
    if (ne == 0) return;

    // Find highest core and lowest valence energies
    int ic = 0;
    int iv = ne + 1;
    for (int ie = 1; ie <= ne; ++ie) {
        if (icv[ie - 1] == -1) {
            ic = ie;
        } else {
            if (ie < iv) iv = ie;
        }
    }

    // Change assignment from core to valence if core state above lowest valence
    for (int ie = iv; ie < ic; ++ie) {
        if (icv[ie] < 0) {
            int iph_local = ip[ie];
            icv[ie] = 1;
            ival[iph_local * (lx + 1) + ll[ie]] = 1;
            // Update occupation number
            xnvmu[iph_local * (lx + 1) + ll[ie]] += 4 * ll[ie] + 2;
            // Update valence density
            int iorb = iiorb[iph_local * (lx + 1) + ll[ie]];
            for (int ir = 0; ir < s251; ++ir) {
                // dgc is [251][30][nphx+2]; access dgc(ir+1, iorb+1, iph)
                int dgc_idx = iph_local * (s251 * 30) + iorb * s251 + ir;
                int dpc_idx = dgc_idx;
                double dgc_val = dgc[dgc_idx];
                double dpc_val = dpc[dpc_idx];
                edenvl[iph_local * s251 + ir] += 2.0 * (ll[ie] + 1) *
                    (dgc_val * dgc_val + dpc_val * dpc_val) / (ri05[ir] * ri05[ir]);
                if (ll[ie] != 0) {
                    int dgc_idx2 = iph_local * (s251 * 30) + (iorb - 1) * s251 + ir;
                    double dgc_val2 = dgc[dgc_idx2];
                    double dpc_val2 = dpc[dgc_idx2];
                    edenvl[iph_local * s251 + ir] += 2.0 * ll[ie] *
                        (dgc_val2 * dgc_val2 + dpc_val2 * dpc_val2) / (ri05[ir] * ri05[ir]);
                }
            }
        }
    }
    ic = iv - 1;

    // Check if suggested ecv is between core and valence
    bool ok = false;
    if (ic > 0) {
        if (iv <= ne) {
            if (ecv - en[ic - 1] > tol && en[iv - 1] - ecv > tol) ok = true;
        } else {
            if (ecv - en[ic - 1] > tol) ok = true;
        }
    } else {
        if (iv <= ne) {
            if (en[iv - 1] - ecv > tol) ok = true;
        }
    }
    if (ok) return;

    // Reassign core states to valence as needed
    while (true) {
        ecv = vint - tol;
        if (iv <= ne) ecv = std::min(ecv, en[iv - 1] - tol);
        if (ic == 0) break;
        if (ecv - en[ic - 1] > tol) break;

        // Reassign last core state to valence
        ic--;
        iv--;
        icv[iv - 1] = 1;
        ival[ip[iv - 1] * (lx + 1) + ll[iv - 1]] = 1;
        xnvmu[ip[iv - 1] * (lx + 1) + ll[iv - 1]] += 4 * ll[iv - 1] + 2;

        // Update valence density
        int iph_local = ip[iv - 1];
        int iorb = iiorb[iph_local * (lx + 1) + ll[iv - 1]];
        for (int ir = 0; ir < s251; ++ir) {
            int dgc_idx = iph_local * (s251 * 30) + iorb * s251 + ir;
            double dgc_val = dgc[dgc_idx];
            double dpc_val = dpc[dgc_idx];
            edenvl[iph_local * s251 + ir] += 2.0 * (ll[iv - 1] + 1) *
                (dgc_val * dgc_val + dpc_val * dpc_val) / (ri05[ir] * ri05[ir]);
            if (ll[iv - 1] != 0) {
                int dgc_idx2 = iph_local * (s251 * 30) + (iorb - 1) * s251 + ir;
                double dgc_val2 = dgc[dgc_idx2];
                double dpc_val2 = dpc[dgc_idx2];
                edenvl[iph_local * s251 + ir] += 2.0 * ll[iv - 1] *
                    (dgc_val2 * dgc_val2 + dpc_val2 * dpc_val2) / (ri05[ir] * ri05[ir]);
            }
        }
    }

    // Update core-valence separation in xnval array
    for (int ie = iv - 1; ie < ne; ++ie) {
        int iph_local = ip[ie];
        int lll = ll[ie];
        int iorb = iiorb[iph_local * (lx + 1) + lll];
        if (xnval[iph_local * 30 + iorb] < 0.1) {
            xnval[iph_local * 30 + iorb] = 2.0 * lll + 2.0;
            if (lll > 0) xnval[iph_local * 30 + (iorb - 1)] = 2.0 * lll;
        }
    }
}

} // namespace feff::pot
