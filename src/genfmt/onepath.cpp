// Single-path F-matrix calculation.
// Converted from GENFMT/onepath.f (~647 lines)
//
// Computes F-matrix for a single scattering path.

#include "onepath.hpp"
#include "genfmt_prep.hpp"
#include "pathgeom.hpp"
#include "setlam.hpp"
#include "rot3i.hpp"
#include "sclmz.hpp"
#include "fmtrxi.hpp"
#include "mmtr.hpp"
#include "mmtrxi.hpp"
#include "xstar.hpp"

#include "../common/logging.hpp"
#include "../common/physics_utils.hpp"
#include "../rdinp/mkptz.hpp"
#include <feff/constants.hpp>
#include <feff/dimensions.hpp>

#include <cmath>
#include <complex>
#include <sstream>
#include <algorithm>

namespace feff::genfmt {

void onepath(const std::string& phpad, int index, int nleg, double deg,
             int iorder, int ipot[], double rat[][legtot + 2],
             int iz_in[],
             int ipol, const double evec[3], double elpty,
             const double xivec[3],
             bool write_feffdat, bool write_xdi, bool verbose,
             OnePathResult& result, OnePathPotInfo& potinfo) {

    constexpr double eps = 1.0e-16;

    // Data structures
    PhaseData pd;
    SpinPhaseData spd;
    NlmData nlm_data;
    RotationMatrixData rm;
    ClmzData clmz;
    FmatrixData fmat;
    LambdaData lam;
    BmatrixData bmatd;
    GenfmtPrepResult prep;

    // Local arrays
    FeffComplex rho[legtot];
    FeffComplex pmati[lamtot][lamtot][2];
    FeffComplex cchi[nex];
    FeffComplex rkk[nex][8];
    FeffComplex ck_local[nex];
    int lind[8];

    // Initialize iz
    for (int i = 0; i <= nphx; ++i) pd.iz[i] = 0;
    for (int i = 0; i <= nphx; ++i) {
        if (iz_in) pd.iz[i] = iz_in[i];
    }

    // Copy ipot into pd
    for (int i = 0; i <= legtot; ++i) {
        pd.ipot[i] = ipot[i];
    }

    // Copy rat into pd
    for (int j = 0; j < 3; ++j)
        for (int i = 0; i <= legtot + 1; ++i)
            pd.rat[j][i] = rat[j][i];

    // Set various flags
    int ispin = 0;
    int le2 = 0;
    double angks = 0.0;
    FeffComplex ptz[3][3] = {};

    // ---- Initialization ----
    genfmt_prep(phpad, ispin, pd, spd, nlm_data, prep);

    potinfo.xmu = prep.xmu;
    potinfo.edge = prep.edge;
    potinfo.rnrmav = prep.rnrmav;
    potinfo.rs = prep.rs;
    potinfo.vint = prep.vint;
    potinfo.xkf = std::abs(prep.ck[prep.ik0]);

    int ne = prep.ne;
    int ne1 = prep.ne1;
    int npot = prep.npot;
    int lmaxp1 = prep.lmaxp1;
    int ilinit = prep.ilinit;
    int kinit = prep.kinit;
    int nsp = prep.nsp;

    // Set gamach
    common::setgam(pd.iz[0], prep.ihole, potinfo.gamach);

    // Copy nleg, set nsc
    pd.nleg = nleg;
    int nsc;

    // Compute path geometry (ri, beta, eta)
    pathgeom(nleg, nsc, ipol, pd.rat, pd.ipot, pd.ri, pd.beta, pd.eta);
    pd.nsc = nsc;

    // Initialize polarization tensor via mkptz (matching Fortran onepath.f line 261)
    {
        double spvec[3] = {0.0, 0.0, 0.0};
        double atarr[1][3] = {{0.0, 0.0, 0.0}};  // dummy for mkptz
        double evec_local[3] = {evec[0], evec[1], evec[2]};
        double xivec_local[3] = {xivec[0], xivec[1], xivec[2]};
        feff::rdinp::mkptz(ipol, elpty, evec_local, xivec_local, ispin, spvec,
                           1, atarr, angks, le2, ptz);
    }

    // Copy outputs to result
    for (int i = 0; i < legtot; ++i) result.ri[i] = pd.ri[i];
    for (int i = 0; i <= legtot; ++i) result.beta[i] = pd.beta[i];
    for (int i = 0; i <= legtot + 1; ++i) result.eta[i] = pd.eta[i];

    // Compute reff
    double reff = 0.0;
    for (int i = 1; i <= nleg; ++i) reff += pd.ri[i];
    reff /= 2.0;

    // Set lambda for low k
    setlam(iorder, 1, pd.beta, nsc, nleg, ilinit, lam);

    // Calculate rotation matrices
    for (int isc = 1; isc <= nleg; ++isc) {
        rot3i(lmaxp1, lam.mmaxp1, isc, pd.beta, rm, isc - 1);
    }
    if (ipol > 0) {
        rot3i(ilinit + 1, ilinit + 1, nleg + 1, pd.beta, rm, nleg);
    }

    // Zero cchi
    for (int ie = 0; ie < ne; ++ie) {
        cchi[ie] = FeffComplex(0.0, 0.0);
    }

    // ---- Spin loop ----
    for (int is = 0; is < nsp; ++is) {
        int ispin_eff = (nsp == 1) ? ispin : is + 1;
        mmtr(bmatd, ipol, ispin_eff, le2, angks, ptz, lind,
             rm, pd.eta, nsc, nleg, kinit, ilinit);

        // Copy spin-specific data
        for (int ie = 0; ie < ne; ++ie) {
            pd.eref[ie] = spd.eref2[ie][is];
        }
        for (int iph = 0; iph <= npot; ++iph) {
            for (int ie = 0; ie < ne; ++ie) {
                for (int il = -pd.lmax[ie][iph]; il <= pd.lmax[ie][iph]; ++il) {
                    ph_at(pd, ie, il, iph) = ph4_at(spd, ie, il, is, iph);
                }
            }
        }
        for (int ie = 0; ie < ne; ++ie) {
            for (int kdif = 0; kdif < 8; ++kdif) {
                rkk[ie][kdif] = spd.rkk2[ie][kdif][is];
            }
        }
        for (int ie = 0; ie < ne; ++ie) {
            ck_local[ie] = std::sqrt(2.0 * (pd.em[ie] - pd.eref[ie]));
        }

        // ---- Energy loop ----
        for (int ie = 0; ie < ne; ++ie) {
            for (int ileg = 1; ileg <= nleg; ++ileg) {
                rho[ileg - 1] = ck_local[ie] * pd.ri[ileg];
            }

            if (std::abs(ck_local[ie]) <= eps) {
                std::ostringstream ss;
                ss << " genfmt: ck=0.  ie, ck(ie) " << ie
                   << " " << ck_local[ie].real() << " " << ck_local[ie].imag();
                common::logger().wlog(ss.str());
                continue;
            }

            // Zero clmi
            for (int il = 0; il < legtot; ++il)
                for (int im = 0; im < mtot + ntot + 1; ++im)
                    for (int il2 = 0; il2 < ltot + 1; ++il2)
                        clmz.clmi[il2][im][il] = FeffComplex(0.0, 0.0);

            int mnmxp1 = lam.mmaxp1 + lam.nmax;

            for (int ileg = 1; ileg <= nleg; ++ileg) {
                int isc0 = ileg - 1;
                if (isc0 == 0) isc0 = nleg;
                int isc1 = ileg;
                int lxp1 = std::max(pd.lmax[ie][pd.ipot[isc0]] + 1,
                                    pd.lmax[ie][pd.ipot[isc1]] + 1);
                int mnp1 = std::min(lxp1, mnmxp1);
                sclmz(rho, lxp1, mnp1, ileg - 1, clmz);
            }

            // First matrix
            fmtrxi(lam.lamx, lam.laml0x, ie, 1, 0,
                   clmz, lam, nlm_data, rm, pd, fmat);

            if (nleg > 2) {
                fmtrxi(lam.laml0x, lam.lamx, ie, nleg - 1, nleg - 2,
                       clmz, lam, nlm_data, rm, pd, fmat);
            }

            for (int ilegp = 2; ilegp <= nsc - 1; ++ilegp) {
                int ileg_idx = ilegp;
                fmtrxi(lam.lamx, lam.lamx, ie, ileg_idx, ilegp - 1,
                       clmz, lam, nlm_data, rm, pd, fmat);
            }

            // Matrix multiplication
            int indp = 0;
            for (int lmp = 0; lmp < lam.laml0x; ++lmp) {
                for (int lm = 0; lm < lam.lamx; ++lm) {
                    pmati[lm][lmp][indp] = fmat.fmati[lm][lmp][0];
                }
            }

            for (int isc = 2; isc <= nleg - 1; ++isc) {
                indp = 1 - (isc % 2);
                int indp0 = 1 - indp;

                for (int lmp = 0; lmp < lam.laml0x; ++lmp) {
                    for (int lm = 0; lm < lam.lamx; ++lm) {
                        FeffComplex pllp(0.0, 0.0);
                        for (int lmi = 0; lmi < lam.lamx; ++lmi) {
                            pllp += fmat.fmati[lm][lmi][isc - 1] *
                                    pmati[lmi][lmp][indp0];
                        }
                        pmati[lm][lmp][indp] = pllp;
                    }
                }
            }

            FeffComplex srho(0.0, 0.0);
            FeffComplex prho(1.0, 0.0);
            for (int ileg = 0; ileg < nleg; ++ileg) {
                srho += rho[ileg];
                prho *= rho[ileg];
            }

            // Termination matrix
            mmtrxi(rkk, lam.laml0x, bmatd, ie, 0, nleg - 1, lind,
                   clmz, lam, nlm_data, pd.eta, fmat);

            // Final trace
            FeffComplex ptrac(0.0, 0.0);
            for (int lm = 0; lm < lam.laml0x; ++lm) {
                for (int lmp = 0; lmp < lam.laml0x; ++lmp) {
                    ptrac += fmat.fmati[lm][lmp][nleg - 1] *
                             pmati[lmp][lm][indp];
                }
            }

            FeffComplex cfac = std::exp(coni * (srho - 2.0 * prep.xk[ie] * reff)) / prho;
            if (nsp == 2 && is == 0) cfac = -cfac;
            cchi[ie] += ptrac * cfac;
        }
        // end energy loop
    }
    // end spin loop

    // Compute phase and amplitude for output
    float phff[nex], amff_arr[nex];
    double phffo = 0.0;
    for (int ie = 0; ie < ne1; ++ie) {
        phff[ie] = 0.0f;
        if (std::abs(cchi[ie]) >= eps) {
            phff[ie] = static_cast<float>(std::atan2(cchi[ie].imag(), cchi[ie].real()));
        }
        if (ie > 0) {
            double ph_d = static_cast<double>(phff[ie]);
            common::pijump(ph_d, phffo);
            phff[ie] = static_cast<float>(ph_d);
        }
        phffo = static_cast<double>(phff[ie]);
        amff_arr[ie] = static_cast<float>(std::abs(cchi[ie]));
    }

    // Store result columns
    // The fdtarr function would compute these columns; we compute them directly.
    result.ne1 = ne1;
    for (int ie = 0; ie < ne1; ++ie) {
        result.col1[ie] = prep.xk[ie];                         // k-grid
        result.col2[ie] = ph_at(pd, ie, prep.ll, 0).real();    // central atom phase
        result.col3[ie] = static_cast<double>(amff_arr[ie]);    // |F_eff|
        result.col4[ie] = static_cast<double>(phff[ie]);        // phase of F_eff
        result.col5[ie] = 1.0;                                  // reduction factor (placeholder)
        result.col6[ie] = (prep.ck[ie].imag() != 0.0) ?
                          1.0 / (2.0 * prep.ck[ie].imag()) : 0.0;  // mean free path
        result.col7[ie] = prep.ck[ie].real();                   // Re(p)
    }
}

} // namespace feff::genfmt
